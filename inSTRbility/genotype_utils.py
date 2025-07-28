from phasing_utils import haplocluster_reads, qvalue_phasing
import numpy as np
from sklearn.cluster import KMeans
import warnings
from threadpoolctl import threadpool_limits


def length_genotyper(locus_key, hallele_counter, LOCI_VARS, read_indices, haploid):
    """
    Build genotype based on the allele lengths of the reads
    Args:
        chrom (str): Chromosome name.
        locus_start (int): Start position of the locus.
        locus_end (int): End position of the locus.
        locus_key (str): Key for the locus.
        hallele_counter (dict): Dictionary with allele lengths and their counts.
        LOCI_INFO (dict): Global locus information.
        LOCI_VARS (dict): Global locus variations.
        read_indices (list): List of read indices.
        ref (pysam.FastaFile): Reference genome.
        out (file): Output file handle.
        haploid (bool): Whether the analysis is for haploid genomes.
    Returns:
        list: [state, fail_code] where state is True if genotype was successfully built,
              and fail_code is the point to skip in the analysis.
    """

    read_indices = sorted(read_indices)

    singleton_alens = [item[0] for item in hallele_counter.items() if item[1]==1] # allele with 1 read contribution
    goodcov_alens   = set(hallele_counter.keys()) - set(singleton_alens) # allele with more than 1 read contribution

    allele_lengths = []
    for rindex in read_indices:
        read_alen = LOCI_VARS[locus_key]['A'][rindex][0]
        if read_alen in singleton_alens: # checking if the '1 read - allele' is nearby any of other 'good read - allele'
            for galen in goodcov_alens:
                if read_alen in range(galen-10, galen+10): # '1 read - allele' is considered if other allele are within 10 bp on either of the side
                    allele_lengths.append(galen)
                    break
        else:
            allele_lengths.append(read_alen)

    if allele_lengths == []:
        return [False, 6]

    # Using KMeans clustering to find two clusters of allele lengths
    data = np.array(allele_lengths)
    data = data.reshape(-1, 1)
    with threadpool_limits(limits=1):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            kmeans = KMeans(n_clusters=2, init='k-means++', n_init=5, random_state=0).fit(data)
    cluster_labels = kmeans.labels_
    C1 = [i for i, x in enumerate(cluster_labels) if x == 0]
    C2 = [i for i, x in enumerate(cluster_labels) if x == 1]

    haplotypes = ([read_indices[idx] for idx in C1], [read_indices[idx] for idx in C2])
    cutoff     = 0.15*len(allele_lengths) # 15%

    # if the chromosome is mentioned as haploid
    if haploid:
        # choosing the homozygous allele length
        haplotypes = (haplotypes[0] + haplotypes[1], [])

    # if chromosome is not haploid and both the clusters have sufficient reads
    elif (C1!=[] and len(C1)>=cutoff) and (C2!=[] and len(C2)>=cutoff):
        reads_phased = 'length-based'
        pass

    # if only the first cluster has sufficient reads
    elif C1!=[] and len(C1)>=cutoff:
        haplotypes = (haplotypes[0] + haplotypes[1], [])
        reads_phased = 'unphased'

    # if only the second cluster has sufficient reads
    elif C2!=[] and len(C2)>=cutoff:
        haplotypes = (haplotypes[0] + haplotypes[1], [])
        reads_phased = 'unphased'

    else:
        reads_phased = 'unphased'; fail_code = 6 # write allele distribution with only one read supporting to it in vcf
        return [reads_phased, fail_code] # write allele distribution with only one read supporting to it in vcf

    haplotags = []
    for read_idx in LOCI_VARS[locus_key]['RDS']:
        if read_idx in haplotypes[0]:   haplotags.append(1)
        elif read_idx in haplotypes[1]: haplotags.append(2)
        else: haplotags.append(None)

    LOCI_VARS[locus_key]['HP'] = haplotags

    fail_code = 10  # default fail code for the locus
    return [reads_phased, fail_code]


def snp_filtering(snp_position, SNPS, args, locus_start, locus_end):
    """
    Filters SNPs based on coverage and distance from the locus.

    Args:
        snp_position (int): The position of the SNP.
        SNPS (dict): Dictionary containing SNP information.
        args (Namespace): Command line arguments.
        locus_start (int): Start position of the locus.
        locus_end (int): End position of the locus.

    Returns:
        bool: True if the SNP passes the filters, False otherwise.
    """

    alt_cov = []
    for allele in SNPS[snp_position]:
        if allele not in ['COV', 'Q', 'r']:
            alt_cov.append(len(SNPS[snp_position][allele]))
    alt_cov = max(alt_cov) if alt_cov else 0

    return (snp_position in SNPS) and (SNPS[snp_position]['COV'] >= 3) and \
           alt_cov > 2 and \
           (locus_start - args.snp_dist < snp_position < locus_end + args.snp_dist)


def assign_haplogroups(chrom, locus_start, locus_end, locus_key, LOCI_INFO, LOCI_VARS, READ_VARS, READ_NAMES,
                       SNPS, hallele_counter, ref, out, args, read_max_limit, haploid, homozygous=False):
    """
    Analyse the genotype for a given locus
    """

    reads_phased = False   # stores the method used for phasing the reads

    read_indices = LOCI_VARS[locus_key]['RDS']
    if read_max_limit: read_indices = read_indices[:args.max_reads]

    if haploid:
        reads_phased, fail_code = length_genotyper(locus_key, hallele_counter, LOCI_VARS, read_indices, haploid)
        return [reads_phased, fail_code]

    # snp positions that are recorded in all the reads near this pocus
    snp_positions = set()
    for rindex in read_indices:
        snp_positions |= (READ_VARS[rindex]['X'])

    # considering only those SNPs that are within the locus range and have sufficient coverage of 3 reads
    snp_positions = sorted(list(filter(lambda x: snp_filtering(x, SNPS, args, locus_start, locus_end), snp_positions)))

    read_indices    = set(read_indices)
    snp_allelereads = {}    # stores the SNP information for each position {pos: {'cov': 0, 'reads': set(), 'alleles': { nuc: {reads} }, 'Q': {}}}
    alt_snp_cov = {}    # coverage for non-ref nucleotides in all SNP positions

    # checking for coverage and average quality for non-reference nucleotides in reads
    for pos in snp_positions:
        alt_alleles   = [allele for allele in SNPS[pos] if allele not in ['COV', 'Q', 'r']]
        allele_coverages = set()
        for allele in alt_alleles:
            allele_reads = SNPS[pos][allele].intersection(read_indices)

            # if the non-reference nucleotide allele is supported only by one read
            if len(allele_reads) <= 1: continue

            average_Q = (sum([SNPS[pos]['Q'][read_idx] for read_idx in allele_reads])/len(allele_reads))

            # the average quality of one nuc is less than 13 then we drop the SNP position
            if average_Q <= 13: # phredQ 13 is equivalent to 0.05 error rate
                continue

            allele_coverages.add(len(allele_reads))

        # only consider the SNP position if the coverage for the nucleotide is sufficient and
        # base quality of the non reference nucleotides is good
        if len(allele_coverages) == 0: continue
        else: alt_snp_cov[pos] = max(allele_coverages)

        snp_allelereads[pos] = { 'cov': 0, 'reads': set(), 'alleles': {}, 'Q': {} }
        for allele in SNPS[pos]:
            if (allele == 'COV') or (allele == 'Q'): continue
            snp_allelereads[pos]['alleles'][allele] = SNPS[pos][allele].intersection(read_indices)
            snp_allelereads[pos]['cov'] += len(snp_allelereads[pos]['alleles'][allele])
            if allele != 'r':
                snp_allelereads[pos]['Q'].update(dict([(read_idx, SNPS[pos]['Q'][read_idx]) for read_idx in snp_allelereads[pos]['alleles'][allele]]))

    reads_phased = 'snp-based'
    haplotypes = set()

    # remove the SNP positions that have less than 5 reads supporting it
    del_positions = list(filter(lambda x: (snp_allelereads[x]['cov'] < 5) or (len(snp_allelereads[x]['alleles']) == 1), snp_allelereads.keys()))
    for pos in del_positions: del snp_allelereads[pos]

    if len(snp_allelereads) == 0: # if there are no SNPs left after filtering
        fail_code = 2
        reads_phased = 'unphased'
        return [reads_phased, fail_code]

    else:
        # ordering the SNPS based on coverage
        coverage_sorted_snps = sorted(snp_allelereads.keys(), key = lambda item : alt_snp_cov[item], reverse = True)
        haplotypes, fail_code, num_SNPs = haplocluster_reads(coverage_sorted_snps, snp_allelereads, read_indices, SNPS, READ_NAMES, args) # SNP ifo and supporting reads for specific locus are given to the phasing function

    if fail_code == 2 or fail_code == 3: # if the SNPs are not sufficient to phase the reads
        reads_phased = 'unphased'

    if haplotypes == (): # if the loci has no significant snps
        if homozygous:
            reads_phased = 'unphased'
        else:
            reads_phased, fail_code = length_genotyper(locus_key, hallele_counter, LOCI_VARS, read_indices, haploid)

        return [reads_phased, fail_code]

    haplotags = []
    for read_idx in LOCI_VARS[locus_key]['RDS']:
        if read_idx in haplotypes[0]:   haplotags.append(1)
        elif read_idx in haplotypes[1]: haplotags.append(2)
        else: haplotags.append(None)

    LOCI_VARS[locus_key]['HP'] = haplotags

    return [reads_phased, fail_code]
