# from locus_utils import process_locus
from cstag_utils import parse_cstag
from cigar_utils import parse_cigar
from operation_utils import update_homopolymer_coords
from locus_utils import locus_processor

import sys
import pysam
import numpy as np


def count_loci(tbx, chrom, first_region, last_region):
    """
    Counts the number of loci in a given chromosome range.
    Args:
        tbx (pysam.TabixFile): Tabix file for fetching genomic regions.
        chrom (str): Chromosome name.
        first_region (tuple): Start and end coordinates of the first region.
        last_region (tuple): Start and end coordinates of the last region.
    Returns:
        int: Total number of loci in the specified range.
    """

    total_loci = 0
    for row in tbx.fetch(chrom, first_region[0], last_region[1]):
        row = row.split('\t')
        if (total_loci == 0) and (int(row[1]) != first_region[0] or int(row[2]) != first_region[1]):
            continue
        total_loci += 1
        if last_region[0] == int(row[1]) and last_region[1] == int(row[2]):
            break
    return total_loci


def remove_read(READ_ENDS, READ_IDXS, READ_VARS, SNPS, SORT_SNPS, PREV_READS, READ_NAMES):
    """
    Removes a read from the tracking lists and updates the SNPs and previous reads.
    Args:
        READ_ENDS (list): List of end positions for each read.
        READ_IDXS (list): List of indices for each read.
        READ_VARS (dict): Dictionary containing variations for each read.
        SNPS (dict): Dictionary containing SNP positions and their coverage.
        SORT_SNPS (list): Sorted list of SNP positions.
        PREV_READS (set): Set of previously processed reads.
    """

    popped      = READ_ENDS.pop(0)
    popped_ridx = READ_IDXS.pop(0)
    if popped_ridx in READ_VARS:
        # remove the SNP information associated with the read
        for pos in READ_VARS[popped_ridx]['X']:
            if pos in SNPS: SNPS[pos]['COV'] -= 1
        del_snps = []
        for pos in SNPS: #!! modify this to terminate iteration after some condition
            if SNPS[pos]['COV'] == 0:
                del_snps.append(pos)
                SORT_SNPS.remove(pos)

        for snp in del_snps: del SNPS[snp]
        del READ_VARS[popped_ridx]; del del_snps
        del READ_NAMES[popped_ridx]

        if popped_ridx in PREV_READS: PREV_READS.remove(popped_ridx)


def remove_locus(LOCUS_ENDS, LOCUS_KEYS, LOCI_INFO, LOCI_VARS):
    """
    Removes a read from the tracking lists and updates the SNPs and previous reads.
    Args:
        READ_ENDS (list): List of end positions for each read.
        READ_IDXS (list): List of indices for each read.
        READ_VARS (dict): Dictionary containing variations for each read.
        SNPS (dict): Dictionary containing SNP positions and their coverage.
        SORT_SNPS (list): Sorted list of SNP positions.
        PREV_READS (set): Set of previously processed reads.
    """

    popped      = LOCUS_ENDS.pop(0)
    locus_key = LOCUS_KEYS.pop(0)
    del LOCI_VARS[locus_key]
    del LOCI_INFO[locus_key]


def remove_snp(SORT_SNPS, SNPS, READ_VARS):
    """
    Removes a SNP from the sorted list and updates the SNP dictionary and read variations.
    Args:
        SORT_SNPS (list): Sorted list of SNP positions.
        SNPS (dict): Dictionary containing SNP positions and their coverage.
        READ_VARS (dict): Dictionary containing variations for each read.
    """

    first_SNP = SORT_SNPS[0]
    first_SNP_reads = set()
    for nuc in SNPS[first_SNP]:
        if nuc == 'COV' or nuc == 'Q': continue
        first_SNP_reads |= SNPS[first_SNP][nuc]
    drop_SNP = True
    for _ in first_SNP_reads:
        if _ in READ_VARS:
            drop_SNP = False; break

    if drop_SNP:
        # if the SNP is not useful for the current read then we remove it
        del SNPS[first_SNP]
        SORT_SNPS.pop(0)
        return True
    else:
        # if the SNP is useful for the current read then we stop removing SNPs
        return False


def record_readloci(args, tbx, read, first_region, last_region, LOCI_INFO, LOCI_ENDS,
                    LOCI_KEYS, LOCI_VARS, loci_coords, loci_keys, read_loci_variations):
    """
    Records the variation of a read at a specific position.
    Args:
        args (argparse.Namespace): Parsed command line arguments.
        tbx (pysam.TabixFile): Tabix file for fetching genomic regions.
        read (pysam.AlignedSegment): The read to be processed.
        first_region (tuple): Start and end coordinates of the first region.
        last_region (tuple): Start and end coordinates of the last region.
        LOCI_INFO (dict): Dictionary containing information about each locus.
        LOCI_ENDS (list): List of end positions for each locus.
        LOCI_KEYS (list): List of locus keys.
        LOCI_VARS (dict): Dictionary containing variations for each locus.
        loci_coords (list): List of coordinates for each locus.
        loci_keys (list): List of keys for each locus.
        read_loci_variations (dict): Dictionary to store variations for the read.
    """
    read_chrom = read.reference_name
    read_start = read.reference_start
    read_end   = read.reference_end

    for row in tbx.fetch(read_chrom, read_start, read_end):

        # adjust read start and end based on soft and hard clippings
        # soft and hard clippings do not consume the reference bases

        row = row.split('\t')
        locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start

        # if the locus is within the range of coordinates assigned to the thread
        if (first_region[0] <= locus_start) and (locus_end <= last_region[1]):
            if locus_start == first_region[0]:
                if locus_end == first_region[1]: pass
                # checking if the locus actually matches the first region
                else: continue
            pass
        elif locus_start < first_region[0]:
            # if the locus is before the first region
            continue
        elif locus_start >= last_region:
            # if the locus is beyond the last region
            break

        # if only the read completely covers the repeat
        if (read_start <= locus_start) and (locus_end <= read_end):

            # the flank distances on the left and right of the repeat to be considered
            left_flank = args.flank; right_flank = args.flank

            # adjust the flank distances based on the read start
            if (locus_start - args.flank) < read_start:
                left_flank = locus_start - read_start
            if (locus_end + args.flank) > read_end:
                right_flank  = read_end - locus_end

            loci_coords.append((locus_start, locus_end, left_flank, right_flank))

            locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
            loci_keys.append(locus_key)
            read_loci_variations[locus_key] = {'HALEN': locus_len, 'ALEN': locus_len, 'RLEN': locus_len, 'S': [],
                                               'tracked': False, 'range': [-1, -1]}

            if locus_key not in LOCI_VARS:
                LOCI_VARS[locus_key] = {'RLEN': locus_len, 'RDS': [], 'A': {}, 'S': {}, 'RNG': {}, 'HP': [], 'Q': []}
                LOCI_INFO[locus_key] = row

                # adding the locus key when it is first encountered
                LOCI_ENDS.append(locus_end)
                LOCI_KEYS.append(locus_key)


def extract_reads(args, bam_file, coordinate_range, tidx, karyotype):
    """
    The first function of the inSTRbility pipeline. For a set of coordinate ranges
    this fcuntions processes all the reads and genotypes all the loci in the coordinate range.
    Args:
        args (argparse.Namespace): Parsed command line arguments.
        bam_file (str): Path to the BAM file containing aligned reads.
        coordinate_range (list): List of tuples containing chromosome, start, and end coordinates for each contig.
        tidx (int): Thread index for parallel processing.
        karyotype (bool): Indicates if the sample is haploid (e.g., X or Y chromosome).
    Returns:
        None
    """

    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(args.bed)
    ref  = pysam.FastaFile(args.ref)
    bam  = pysam.AlignmentFile(bam_file, args.aln_format)

    # set output file and log file names
    if tidx != -1:
        if tidx==0:
            out = open(f'{args.output}', 'w')
        else:
            out = open(f'{args.output}_thread_{tidx}', 'w')

    else:
        out = open(f'{args.output}', 'w')

    # used for identifying the flanking region within the read incase there's indels nearby
    flank_length = 50

    print('locus_id', 'fail_code', 'zygosity', 'phased_on', 'read_index', 'read_start_pos', 'read_end_pos',
          'allele_length', 'sequence', 'avg_basequal', 'haplotag', sep='\t', file=out)

    # contig range is a list of tuples
    # each tuple containes a chromosome, the coordinates of the first region and the coordinates of the last region
    for region in coordinate_range:

        prev_locus_end = 0

        # the first and last repeat regions in the contig
        chrom, first_region, last_region = region

        haploid = False
        if (chrom in {'chrX', 'chrY', 'X', 'Y'}) and karyotype: haploid = True

        total_loci = count_loci(tbx, chrom, first_region, last_region)
        print(f"> {chrom} {first_region} {last_region} Total loci = ", total_loci, file=sys.stderr)

        # {pos: {'COV': coverage, Q: {read_index: base_phredQ}, 'nuc': {read_idx1, read_idx2}}
        SNPS = dict()           # tracking the encountered SNPs
        # {read_index: {'S': start, 'E': end, 'X': set(), 'D': set()}}
        READ_VARS = dict()      # tracking the variations on each read
        # {locus_key: {'RLEN': repeat length, 'RDS': [read_index1, read_index2], 'A': {read_index: [halen, alen]},
        #  'S': {read_index: seq}, 'RNG': {read_index: [start, end]}, 'HP': [haplotag1, haplotag2]}}
        LOCI_VARS = dict()      # tracking the variation for each locus
        # {locus_key: [chrom, start, end, motif]}
        LOCI_INFO = dict()      # saving the information of each loci

        # tracking the loci
        # locus ends, keys | read ends and indices stored in the order of encountering
        LOCI_ENDS = []; LOCI_KEYS = []
        READ_ENDS = []; READ_IDXS = []
        SORT_SNPS = []
        READ_NAMES = dict()

        # set of the previous reads that are processed for the last locus
        PREV_READS = set()

        start_coord = first_region[0]; end_coord = last_region[1]

        # processing the reads within the coordinate range
        read_index = 0
        for read in bam.fetch(chrom, start_coord, end_coord):

            # skip read with low mapping quality
            if read.mapping_quality < args.mapq:
                continue

            read_chrom = read.reference_name
            read_start = read.reference_start
            read_end   = read.reference_end

            while len(LOCI_ENDS) > 0 and read_start > LOCI_ENDS[0]:
                # if the current read is beyond the start of the first locus
                # we process the first locus for its variations
                prev_locus_end = locus_processor(args, tbx, ref, out, LOCI_KEYS, LOCI_ENDS, LOCI_VARS, LOCI_INFO, READ_VARS,
                                                 PREV_READS, READ_NAMES, SNPS, SORT_SNPS, haploid, prev_locus_end, flank_length)

            while len(LOCI_ENDS) > 0 and read_start > READ_ENDS[0]:
                # if the current read is beyond the end of the tracked read
                # we check if the read is useful for processsing any loci and get rid of it
                if len(LOCI_ENDS) > 0 and READ_ENDS[0] > LOCI_ENDS[0]:
                    # if the initial read useful for the first locus being tracked then it is retained
                    break

                else:
                    # remove the read information if the current read is beyond the first read and the locus
                    remove_read(READ_ENDS, READ_IDXS, READ_VARS, SNPS, SORT_SNPS, PREV_READS, READ_NAMES)

            while len(SORT_SNPS) > 0 and read_start > SORT_SNPS[0]:
                # if the current read is beyond the first SNP position
                # we remove the SNP information that is not useful for the current read
                if not remove_snp(SORT_SNPS, SNPS, READ_VARS):
                    # if the SNP is useful for the current read then we stop removing SNPs
                    break

            # if the read is beyond the last locus in the bed file the loop stops
            if read_start > end_coord:
                while len(LOCI_ENDS) > 0:
                    del LOCI_ENDS[0]
                # process the loci left in global_loci_variation
                break

            # set of homopolymer positions within the reference part that is covered by the read
            homopoly_positions = {}

            # repeat loci covered by the read
            loci_coords = []; loci_keys = []
            # information locally saved for each read and all the loci it covers
            read_loci_info = {}
            record_readloci(args, tbx, read, first_region, last_region, LOCI_INFO, LOCI_ENDS,
                            LOCI_KEYS, LOCI_VARS, loci_coords, loci_keys, read_loci_info)

            # if no repeats are covered by the read
            if len(loci_coords) == 0: continue

            read_index += 1
            READ_NAMES[read_index] = read.query_name

            # getting the haplotag information for the read
            haplotag = None
            if (args.haplotag is not None) and read.has_tag(args.haplotag):
                haplotag = read.get_tag(args.haplotag)

            READ_ENDS.append(read_end)
            READ_IDXS.append(read_index)
            READ_VARS[read_index] = {'S': read_start, 'E': read_end, 'X': set(), 'D': set()}

            # parsing through the read CS tag or CIGAR string to identify variation falling within the repeat loci
            if read.has_tag('cs'):
                cs_tag = read.get_tag('cs')
                parse_cstag(read_index, read, loci_keys, loci_coords, read_loci_info, homopoly_positions,
                            READ_VARS, SNPS, SORT_SNPS, flank_length, haploid)
            else:
                parse_cigar(read_index, read, loci_keys, loci_coords, read_loci_info, homopoly_positions,
                            READ_VARS, SNPS, SORT_SNPS, flank_length, haploid, ref)

            for locus_key in read_loci_info:
                start_idx = read_loci_info[locus_key]['range'][0]
                end_idx   = read_loci_info[locus_key]['range'][1]
                length    = end_idx - start_idx
                LOCI_VARS[locus_key]['RDS'].append(read_index)
                LOCI_VARS[locus_key]['HP'].append(haplotag)
                LOCI_VARS[locus_key]['A'][read_index] = [read_loci_info[locus_key]['HALEN'], read_loci_info[locus_key]['ALEN']]
                LOCI_VARS[locus_key]['S'][read_index] = read_loci_info[locus_key]['S']
                LOCI_VARS[locus_key]['RNG'][read_index] = read_loci_info[locus_key]['range']
                LOCI_VARS[locus_key]['RLEN'] = read_loci_info[locus_key]['RLEN']
                base_qualities = read.query_qualities[start_idx:end_idx] if read.query_qualities else None
                if base_qualities is not None and len(base_qualities) > 0:
                    LOCI_VARS[locus_key]['Q'].append(round(np.mean(base_qualities), 2))
                else:
                    LOCI_VARS[locus_key]['Q'].append(None)

        while len(LOCI_ENDS) > 0:
            # if the current read is beyond the start of the first locus
            # we process the first locus for its variations
            prev_locus_end = locus_processor(args, tbx, ref, out, LOCI_KEYS, LOCI_ENDS, LOCI_VARS, LOCI_INFO, READ_VARS,
                                                PREV_READS, READ_NAMES, SNPS, SORT_SNPS, haploid, prev_locus_end, flank_length)

        while len(LOCI_ENDS) > 0:
            # if the current read is beyond the start of the first locus
            # we process the first locus for its variations
            remove_locus(LOCI_ENDS, LOCI_KEYS, LOCI_INFO, LOCI_VARS) # removing the first locus from the tracking lists


    bam.close()
    ref.close()
    tbx.close()
    out.close()
