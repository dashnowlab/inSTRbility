import bisect
from operation_utils import match_jump, deletion_jump, insertion_jump


def parse_cstag(read_index, read, loci_keys, loci_coords, read_loci_variations,
                homopoly_positions, READ_VARS, SNPS, SORT_SNPS, flank_lengthgth, haploid):
    """
    Parse the CS tag for a read and record the variations observed for the read also for the loci
    """

    # initialising sorted SNP positions
    if SORT_SNPS == None: SORT_SNPS = []

    # CS tag operations set
    operations = {':', '-', '+', '*', '=', '~'}
    
    read_start = read.reference_start
    rpos = read.reference_start   # NOTE: The coordinates are 1 based in SAM; pysam uses 0 based coordinates.
    qpos = 0            # starts from 0 the sub string the read sequence in python

    # stores the each locus information in this read 
    loci_info = {}
    for key in loci_keys:
        loci_info[key] = { 'range': [0,0], 'flank_range': [0, 0], 'tracked': False, 'flank_tracked': [False, False] }

    repeat_index = 0
    cigar_one = read.cigartuples[0]

    # if the first cigar operation is soft clip
    # check if the soft clip is close to start of a locus
    if cigar_one[0] == 4:
        soft_clip_length = cigar_one[1]
        for coord in loci_coords:
            if coord[0] - read_start <= flank_lengthgth:
                detect_flank()
                pass
        qpos += cigar_one[1]

    cs_tag = read.get_tag('cs')
    i = 0; cs_len = len(cs_tag)
    while i < cs_len:

        # the sequence match followed by length of the match
        if cs_tag[i] == ':':
            match_len = '0'; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                match_len += cs_tag[i]; i += 1

            match_len = int(match_len)
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(rpos, qpos, match_len, repeat_index, loci_coords, loci_keys, loci_info)

        elif cs_tag[i] == '=':      # sequence match in long CS is followed by nucs which are matching       
            match_len = 0
            while i < cs_len and cs_tag[i] not in operations:
                match_len += 1; i += 1

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(rpos, qpos, match_len, repeat_index, loci_coords, loci_keys, loci_info)

        elif cs_tag[i] == '*':      # substitution of a base; is followed by reference and substituted base
            ref_nuc = cs_tag[i+1]; sub_nuc = cs_tag[i+2]
            i += 3

            if not haploid:
                base_phredQ = read.query_qualities[qpos]
                READ_VARS[read_index]['snps'].add(rpos)
                
                if rpos not in SNPS:
                    SNPS[rpos] = { 'C': 1, sub_nuc: {read_index}, 'Q': { read_index: base_phredQ } }
                    bisect.insort(SORT_SNPS, rpos)
                else:
                    SNPS[rpos]['C'] += 1
                    SNPS[rpos]['Q'][read_index] = base_phredQ
                    if sub_nuc in SNPS[rpos]: 
                        SNPS[rpos][sub_nuc].add(read_index)
                    else: SNPS[rpos][sub_nuc] = {read_index}

            qpos += 1; rpos += 1; match_len = 1
            repeat_index += match_jump(rpos, qpos, match_len, repeat_index, loci_coords, loci_keys, loci_info)

        elif cs_tag[i] == '+':      # insertion; is followed by the inserted bases
            insert = ''; insert_length = 0; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                insert += cs_tag[i]; insert_length += 1
                i += 1
            qpos += insert_length
            repeat_index += insertion_jump(insert_length, insert, rpos, repeat_index, loci_keys, loci_coords,
                                           homopoly_positions, read_loci_variations, qpos)

        elif cs_tag[i] == '-':      # deletion; is followed by the deleted bases
            deletion = ''; deletion_length = 0; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                deletion += cs_tag[i]; deletion_length += 1
                i += 1
            if not haploid:
                READ_VARS[read_index]['dels'] |= set(range(rpos, rpos+deletion_length))
            rpos += deletion_length
            repeat_index += deletion_jump(deletion_length, rpos, repeat_index, loci_keys, loci_coords,
                                          homopoly_positions, read_loci_variations, qpos)
