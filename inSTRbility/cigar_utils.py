import bisect
import sys
import numpy as np
from operation_utils import match_jump, deletion_jump, insertion_jump
from md_utils import parse_mdtag


# !/usr/bin/env python3
def get_substituion_positons(ref, query):
    """
    Returns the indices of the substitution positions between the reference and query sequences.
    Args:
        ref (str): Reference sequence
        query (str): Query sequence 
    Returns:
        list: List of indices where substitutions occur
    """

    if '=' in query:
        array = np.frombuffer(query.encode(), dtype=np.byte)
        substitution_indices = np.where(array != ord('='))[0]
        
    else:
        array1 = np.frombuffer(ref.encode(), dtype=np.byte)
        array2 = np.frombuffer(query.encode(), dtype=np.byte)
        substitution_indices = np.where(array1 != array2)[0]
        
    return substitution_indices.tolist()


def parse_cigar(read_index, read, loci_keys, loci_coords, read_loci_info, homopoly_positions,
                READ_VARS, SNPS, SORT_SNPS, flank_length, haploid, ref):
    """
    Parses the CIGAR tag of a read and updates the read variations and SNPs.
    Args:
        read_index (int): Index of the read in the global read variations
        read (pysam.AlignedSegment): Pysam read object
        loci_keys (list): List of locus keys
        loci_coords (list): List of locus coordinates
        read_loci_info (dict): Dictionary to store read variations for each locus
        homopoly_positions (dict): Dictionary to store homopolymer positions
        READ_VARS (list): List to store global read variations
        SNPS (dict): Dictionary to store global SNP positions
        SORT_SNPS (list): Sorted list of global SNP positions
        flank_length (int): Length of the flanking region to consider
        haploid (bool): Whether the analysis is for haploid data
        ref (pysam.FastaFile): Reference genome file
    Returns:
        None: Updates the global and local read variation and global snp information in place.
    """

    # initialising sorted SNP positions
    if SORT_SNPS == None: SORT_SNPS = []

    chrom = read.reference_name

    rpos = read.reference_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0            # starts from 0 the sub string the read sequence in python

    loci_info = {}
    for key in loci_keys:
        loci_info[key] = { 'range': [0,0], 'flank_range': [0, 0], 'tracked': False, 'flank_tracked': [False, False] }

    repeat_index = 0

    X_tag = False
    insertion_point = {}

    # checking if the read has an MD tag
    md = False
    if read.has_tag('MD'): md = True

    if SORT_SNPS == None: SORT_SNPS = []

    for c, cigar in enumerate(read.cigartuples):
        if cigar[0] == 4:   # soft-clipped read
           qpos += cigar[1]

        elif cigar[0] == 2:     # deletion
            deletion_length = cigar[1]
            if not haploid:
                # recording the deletion positions in the read
                READ_VARS[read_index]['D'] |= set(range(rpos, rpos+deletion_length))
            rpos += cigar[1]
            repeat_index += deletion_jump(rpos, qpos, deletion_length, repeat_index, loci_keys, loci_coords,
                                          read_loci_info, homopoly_positions)

        elif cigar[0] == 1:     # insertion
            insert_length = cigar[1]
            insert = read.query_sequence[qpos:qpos+insert_length]
            insertion_point[rpos] = insert_length
            qpos += insert_length
            repeat_index += insertion_jump(rpos, qpos, insert_length, insert, repeat_index, loci_keys, loci_coords,
                                           read_loci_info, read, ref, homopoly_positions, flank_length, read.reference_start,
                                           read.reference_end, read.cigartuples)

        elif cigar[0] == 0: # match (both equals & difference)
            match_length = cigar[1]
            if not md:
                ref_sequence   = ref.fetch(chrom, rpos, rpos+match_length)
                query_sequence = read.query_sequence[qpos:qpos+match_length]
                substitution_positions = []

                if len(ref_sequence) == len(query_sequence):
                    substitution_positions = get_substituion_positons(ref_sequence, query_sequence)
                else:
                    print('Error in fetching in sequences of ref & read for substitution')
                    sys.exit()

                for sub_pos in substitution_positions:
                    sub_nuc = query_sequence[qpos + sub_pos]
                    base_phredQ = read.query_qualities[qpos + sub_pos]
                    READ_VARS[read_index]['X'].add(rpos + sub_pos)
                    if rpos + sub_pos not in SNPS:
                        SNPS[rpos + sub_pos] = { 'COV': 1, sub_nuc: {read_index}, 'Q': { read_index: base_phredQ } }
                        bisect.insort(SORT_SNPS, rpos + sub_pos)
                    else:
                        SNPS[rpos+sub_pos]['COV'] += 1
                        SNPS[rpos+sub_pos]['Q'][read_index] = base_phredQ
                        if sub_nuc in SNPS[rpos+sub_pos]: 
                            SNPS[rpos + sub_pos][sub_nuc].add(read_index)
                        else: SNPS[rpos + sub_pos][sub_nuc] = {read_index}

            qpos += match_length; rpos += match_length
            repeat_index += match_jump(rpos, qpos, match_length, repeat_index, loci_keys, loci_coords, read_loci_info)

        elif cigar[0] == 7: # exact match (equals)
            X_tag = True
            match_length = cigar[1]; qpos += match_length; rpos += match_length
            repeat_index += match_jump(rpos, qpos, match_length, repeat_index, loci_keys, loci_coords, read_loci_info)

        elif cigar[0] == 8: # substitution (difference)
            X_tag = True
            substitution_length = cigar[1]
            if not haploid:
                for pos in range(substitution_length):
                    sub_nuc = read.query_sequence[qpos + pos]
                    Qval    = read.query_qualities[qpos + pos]
                    READ_VARS[read_index]['X'].add(rpos + pos)
                    if rpos + pos not in SNPS:
                        SNPS[rpos + pos] = { 'COV': 1, sub_nuc: { read_index }, 'Q': { read_index : Qval } }
                        bisect.insort(SORT_SNPS, rpos + pos)
                    else:
                        SNPS[rpos + pos]['COV'] += 1
                        SNPS[rpos + pos]['Q'][read_index] = Qval
                        if sub_nuc in SNPS[rpos + pos]: 
                            SNPS[rpos + pos][sub_nuc].add(read_index)
                        else: SNPS[rpos + pos][sub_nuc] = {read_index}
            qpos += substitution_length; rpos += substitution_length
            repeat_index += match_jump(rpos, qpos, substitution_length, repeat_index, loci_keys, loci_coords, read_loci_info)

    if not X_tag :
        if read.has_tag('MD'):
            if read.cigartuples[0][0] == 4: qpos = read.cigartuples[0][1]
            else: qpos=0
            MD_tag = read.get_tag('MD')
            parse_mdtag(MD_tag, qpos, read.reference_start, READ_VARS, SNPS, read_index,
                        read.query_qualities, read.query_sequence, SORT_SNPS, insertion_point)

    for idx, loc_key in enumerate(loci_keys):
        spos = read_loci_info[loc_key]['range'][0]
        epos = read_loci_info[loc_key]['range'][1]
        read_loci_info[loc_key]['HALEN'] = epos - spos
        read_loci_info[loc_key]['ALEN'] = epos - spos
        read_loci_info[loc_key]['S'] = read.query_sequence[spos:epos]
