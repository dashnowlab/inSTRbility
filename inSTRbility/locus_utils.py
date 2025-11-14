# from realignment_utils import *
from genotype_utils import assign_haplogroups
from collections import Counter
from instabillity_utils import consensus_allele
import sys


def get_mode(list):
    """
    Returns the mode of a list.
    If there are multiple modes, returns the first one found.
    """
    if not list: return None
    return Counter(list).most_common(1)[0][0]


def count_alleles(locus_key, read_indices, LOCI_VARS, allele_counter, hallele_counter):
    """
    Counts the read distribution for each allele length
    """
    for rindex in read_indices:
        halen, alen = LOCI_VARS[locus_key]['A'][rindex]

        try: allele_counter[alen] += 1
        except KeyError: allele_counter[alen] = 1

        try: hallele_counter[halen] += 1
        except KeyError: hallele_counter[halen] = 1


def record_reference_snps(read_indices, locus_start, locus_end, snp_dist, prev_locus_end, READ_VARS, SNPS, SORT_SNPS):
    """
    Records the SNP position as reference for the reads the SNP positions is not present
    Args:
        read_indices (list): List of read indices that cover the locus
        old_reads (set): Set of read indices that have been previously processed
        new_reads (set): Set of read indices that are new at this locus
        READ_VARS (dict): Dictionary containing global read variations
        SNPS (dict): Dictionary containing global SNP positions
        SORT_SNPS (list): Sorted list of SNP positions
    Returns:
        None: Updates the SNPS dictionary with read coverage information
    """

    read_indices = set(read_indices)
    range_start = locus_start - snp_dist
    range_end   = locus_end + snp_dist
    prev_locus_end += snp_dist

    for rindex in read_indices:

        for pos in SORT_SNPS:
            if pos < prev_locus_end: continue
            if not (range_start <= pos <= range_end): continue
            if pos < READ_VARS[rindex]['S']: continue
            if pos > READ_VARS[rindex]['E']: break

            if (pos not in READ_VARS[rindex]['X']) and pos not in READ_VARS[rindex]['D']:
                if 'r' in SNPS[pos]: SNPS[pos]['r'].add(rindex)
                else: SNPS[pos]['r'] = { rindex }
                SNPS[pos]['COV'] += 1


def inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set):
    for locus in near_by_loci:
        if locus[0] <= ins_rpos <= locus[1]:
            sorted_global_ins_rpos_set.add(ins_rpos)
            return 1
    return 0


def process_locus(locus_start, locus_end, locus_key, prev_locus_end, LOCI_VARS,
                  READ_VARS, SNPS, SORT_SNPS, PREV_READS, args):
    """
    processes the locus to phase the reads and get a consensus sequence of the repeat locus.

    Args:
        locus_start (int): Start position of the locus.
        locus_end (int): End position of the locus.
        locus_key (str): Key for the locus in LOCI_VARS.
        prev_locus_end (int): End position of the previous locus.
        LOCI_VARS (dict): Dictionary containing locus variations.
        READ_VARS (dict): Dictionary containing read variations.
        SNPS (dict): Dictionary containing SNP positions and their coverage.
        SORT_SNPS (list): Sorted list of SNP positions.
        PREV_READS (set): Set of previously processed reads.
        args (Namespace): Command line arguments.
    
    Returns:
        tuple: Contains the updated PREV_READS, phasing category, hallele_counter,
               fail_code, read_max_limit, and haplotypes.
    """

    phasing_category, haplotypes = [None, None]

    read_indices = LOCI_VARS[locus_key]['RDS']   # the read indices which cover the locus
    haplotags = LOCI_VARS[locus_key]['HP']
    check = False
    if None in haplotags: check = True
    total_reads = len(read_indices)                             # total number of reads

    read_max_limit = 0

    # remove if the locus has poor coverage
    if total_reads < args.min_reads:
        # coverage of the locus is low
        PREV_READS = set(read_indices)
        fail_code = 1
        return [PREV_READS, phasing_category, {}, fail_code, read_max_limit, haplotypes]

    elif total_reads > args.max_reads:
        # coverage of the locus is high
        # read_indices = read_indices[:args.max_reads]
        # haplotags = haplotags[:args.max_reads]
        read_max_limit = 1

    # reads relevant to this locus
    current_reads = set(read_indices)

    # recording the counts of each allele length across all reads
    allele_counter = {};  hallele_counter = {}
    count_alleles(locus_key, read_indices, LOCI_VARS, allele_counter, hallele_counter)

    hap_status = False
    if args.haplotag is not None: # if the reads are haplotagged we fetch the haplotype of the reads
        haplotypes = ([], [])
        for idx,i in enumerate(haplotags):
            if   i == 1: haplotypes[0].append(read_indices[idx])
            elif i == 2: haplotypes[1].append(read_indices[idx])
        # if both of the haplotypes have reads assigned to them
        hap_status = all([len(hap)>0 for hap in haplotypes])

    # processing haplotagged reads to write into vcf_heterozygous
    # if hap_status & ((haplotags.count(None)/total_reads) <= 0.15):
    if hap_status & (total_reads - haplotags.count(None) >= 20):
        phasing_category = 3 # phased

    elif len(hallele_counter) == 1:
        phasing_category = 1 # homozygous
        homozygous_allele = list(hallele_counter.keys())[0]

    else:
        # removing allele lengths represented by only one read
        filtered_alleles = list(filter(lambda x: hallele_counter[x] > 1, hallele_counter.keys()))
        # if the number of allele lengths left after filtering is only one; assign the locus as homozygous
        if len(filtered_alleles) == 1 and hallele_counter[filtered_alleles[0]]/total_reads >= 0.75:
            phasing_category = 1 # homozygous
            homozygous_allele = filtered_alleles[0]
            reads_of_homozygous = [rindex for rindex in LOCI_VARS[locus_key]['A'] if homozygous_allele == LOCI_VARS[locus_key]['A'][rindex][0]]
        else:
            phasing_category = 2 # ambiguous

    record_reference_snps(read_indices, locus_start, locus_end, args.snp_dist, prev_locus_end, READ_VARS, SNPS, SORT_SNPS)

    PREV_READS = current_reads.copy()
    fail_code = 10  # default fail code for the locus
    return [PREV_READS, phasing_category, hallele_counter, fail_code, read_max_limit, haplotypes]


def locus_processor(args, tbx, ref, out, reads_out, LOCI_KEYS, LOCI_ENDS, LOCI_VARS, LOCI_INFO, READ_VARS,
                    PREV_READS, READ_NAMES, SNPS, SORT_SNPS, haploid, prev_locus_end, flank):
    """
    Processing the information of a locus after processing all the reads that overlap the locus
    Args:
        args (Namespace): Command line arguments
        tbx (pysam.TabixFile): Tabix file object for the reference genome
        ref (pysam.FastaFile): Fasta file object for the reference genome
        LOCI_KEYS (list): List of locus keys
        LOCI_ENDS (list): List of locus end positions
        LOCI_VARS (dict): Dictionary containing locus variations
        LOCI_INFO (dict): Dictionary containing locus information
        READ_VARS (dict): Dictionary containing read variations
        PREV_READS (set): Set of previously processed reads
        SNPS (dict): Dictionary containing SNP positions and their coverage
        SORT_SNPS (list): Sorted list of SNP positions
        haploid (bool): Whether the analysis is haploid or not
        flank (int): Flanking region length around the locus
    Returns:
        prev_locus_end (int): the end position of the last processed locus
    """

    # removes the locus key for from LOCI_KEYS, LOCI_ENDS
    popped    = LOCI_ENDS.pop(0)
    locus_key = LOCI_KEYS.pop(0)
    locus_chrom = LOCI_INFO[locus_key][0]
    locus_start = int(LOCI_INFO[locus_key][1])
    locus_end   = int(LOCI_INFO[locus_key][2])
    motif_length = len(LOCI_INFO[locus_key][3])

    nucleotides = set(['A', 'T', 'C', 'G', 'a', 't', 'c', 'g'])
    if len(set(LOCI_INFO[locus_key][3]) - nucleotides) > 0:
        # if the motif is not a nucleotide sequence, it is a random or complex motif
        # so we set the motif length to -1
        motif_length = -1

    # if LOCI_INFO[locus_key][3] == 'random' or LOCI_INFO[locus_key] == 'complex': motif_length = -1

    near_by_loci = []
    for row in tbx.fetch(locus_chrom, locus_start-flank, locus_end+flank):
        row = row.split('\t')
        near_by_loci.append(( int(row[1]), int(row[2]) ))

    zygosity = ''

    # if the locus has recorded information from the reads
    if locus_key in LOCI_VARS:
        reads_phased = 'unphased'

        if not SORT_SNPS: SORT_SNPS = sorted(SNPS.keys())
        result = process_locus(locus_start, locus_end, locus_key, prev_locus_end, LOCI_VARS,
                               READ_VARS, SNPS, SORT_SNPS, PREV_READS, args)

        PREV_READS, phasing_category, hallele_counter, fail_code, read_max_limit, haplotypes = result

        if phasing_category == 1: # phasing category is homozygous
            zygosity = 'homozygous'
            reads_phased, fail_code = assign_haplogroups(locus_chrom, locus_start, locus_end, locus_key, LOCI_INFO, LOCI_VARS, READ_VARS, READ_NAMES,
                                                         SNPS, hallele_counter, args, read_max_limit, haploid, homozygous=True)

        elif phasing_category == 2:  # phasing category is ambiguous - needs to be phased either based on SNPs of lengths
            reads_phased, fail_code = assign_haplogroups(locus_chrom, locus_start, locus_end, locus_key, LOCI_INFO, LOCI_VARS, READ_VARS, READ_NAMES,
                                                         SNPS, hallele_counter, args, read_max_limit, haploid)

            if reads_phased and reads_phased != 'unphased':
                # based on the haplogroups we identify the allele lengths in the haplogroups
                A = []; B = []
                for _, hp in enumerate(LOCI_VARS[locus_key]['HP']):
                    if hp == 1: A.append(LOCI_VARS[locus_key]['A'][LOCI_VARS[locus_key]['RDS'][_]][0])
                    elif hp == 2: B.append(LOCI_VARS[locus_key]['A'][LOCI_VARS[locus_key]['RDS'][_]][0])
                A = get_mode(A)
                B = get_mode(B)

                # the allele sizes are None if the haplogroup is empty
                if A is None and B is None: zygosity = 'unphased'
                elif A is None or B is None: zygosity = 'homozygous'
                elif A == B: zygosity = 'homozygous'
                else: zygosity = 'heterozygous'
            else:
                fail_messages = {
                    1: 'Locus failed due to less read contribution.',
                    2: 'Locus failed haplotagging due to absence of significant SNPs.',
                    3: 'Locus skipped due to less read contribution in the phased clusters.',
                    6: f'Locus {locus_key} skipped due to wide distribution of alleles with one read supporting to it.',
                }

        elif phasing_category == 3: # phased based on haplotag
            reads_phased = 'haplotag-based'
            A = []; B = []
            for _, hp in enumerate(LOCI_VARS[locus_key]['HP']):
                if hp == 1: A.append(LOCI_VARS[locus_key]['A'][LOCI_VARS[locus_key]['RDS'][_]][0])
                elif hp == 2: B.append(LOCI_VARS[locus_key]['A'][LOCI_VARS[locus_key]['RDS'][_]][0])
            A = get_mode(A)
            B = get_mode(B)

            # the allele sizes are None if the haplogroup is empty
            if A is None and B is None: zygosity = 'unphased'
            elif A is None or B is None: zygosity = 'homozygous'
            elif A == B: zygosity = 'homozygous'
            else: zygosity = 'heterozygous'

        else:
            zygosity = 'unphased'
            pass

        read_indices = LOCI_VARS[locus_key]['RDS']
        sorted_ridx = sorted(list(range(0, len(read_indices))), key=lambda x: LOCI_VARS[locus_key]['HP'][x] if LOCI_VARS[locus_key]['HP'][x] is not None else 3)

        if fail_code == 10 and (reads_phased == 'snp-based' or reads_phased == 'haplotag-based'):
            consensus_allele(locus_key, LOCI_VARS, LOCI_INFO, motif_length, args.pure_stretch, args.num_hpreads, out)

        for ridx in sorted_ridx:
            read_idx = LOCI_VARS[locus_key]['RDS'][ridx]
            fail_messages = {1: 'low-reads', 2: 'nosig-snps', 3: 'low-read-snps', 4: 'low-reads-cluster', 6: 'low-reads', 10: 'pass'}
            if reads_out is not None:
                print(locus_key, fail_messages[fail_code], zygosity, reads_phased, READ_NAMES[read_idx], LOCI_VARS[locus_key]['RNG'][read_idx][0],
                      LOCI_VARS[locus_key]['RNG'][read_idx][1], LOCI_VARS[locus_key]['A'][read_idx][0], LOCI_VARS[locus_key]['S'][read_idx],
                      LOCI_VARS[locus_key]['Q'][ridx], LOCI_VARS[locus_key]['HP'][ridx], sep='\t', file=reads_out)

        del LOCI_VARS[locus_key]
        del LOCI_INFO[locus_key]

    prev_locus_end = popped
    return prev_locus_end
