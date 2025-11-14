#! /usr/bin/env python3
from lib_ssw.pyssw import align_pair


def convert_cigar(cigar):
    """
    Converts the CIGAR string to a list of tuples
    Args:
        cigar: CIGAR string from pysam record

    Returns:
        list of tuples
    """
    cigar_char = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    cigartuples = []
    length = ''
    for c in cigar:
        if c.isdigit(): length += c
        else:
            cigartuples.append((cigar_char[c], int(length)))
            length = ''

    return cigartuples


def convert_cigartuples(cigartuples):
    """
    Converts the CIGAR tuples to a string
    Args:
        cigartuples: CIGAR tuples from pysam record

    Returns:
        CIGAR string
    """
    cigar_char = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}
    cigar = ''
    for c in cigartuples:
        cigar += f'{c[1]}{cigar_char[c[0]]}'

    return cigar


def remove_softclips(cigar):
    """
    Removes the softclips from the CIGAR string

    Args:
        cigar(str): CIGAR string

    Returns:
        CIGAR string without softclips
    """

    new_cigar = ''
    length = ''
    for c in cigar:
        if c.isdigit(): length += c
        else:
            if c != 'S': new_cigar += f'{length}{c}'
            length = ''
    return new_cigar


def complement_cigar(cigar):
    """
    Reverses the CIGAR string
    Args:
        cigar: CIGAR string from pysam record

    Returns:
        list of tuples
    """
    new_cigar = ''
    cigar = remove_softclips(cigar)
    for c in cigar:
        if c == 'I': c = 'D'
        elif c == 'D': c = 'I'
        new_cigar += c
    return new_cigar


def get_indel_position(cigar):
    """
    Get the position of the indel in the CIGAR string
    Args:
        cigar: CIGAR string

    Returns:
        position of the indel in the CIGAR string
    """
    length = ''; pos = 0
    for c in cigar:
        if c.isdigit():
            length += c
        else:
            if c == 'M' or c == '=' or c == 'X':
                pos += int(length)
            elif c == 'I' or c == 'D':
                return [pos, c, int(length)]
            length = ''
    return [-1, -1, -1]


def detect_flank(read, ref, flank_length, locus_start, locus_end, relative_position, indel_type, indel_length, 
                 indel_position, aligned_start, aligned_end, aligned_cigar):
    """
    Identify the flanking sequences of the repeat region in the read
    Args:
        read: pysam read object
        fasta: pysam fasta object
        flank_length: length of the flanking sequence to be considered
        locus_start: start coordinate of the repeat in reference
        locus_end: end coordinate of the repeat in reference

    Returns:
        [start, end, sub_cigar]
        start and end coordinates of the flanking sequence in the read
        sub_cigar: CIGAR string of the repeat alignment to the reference in the read
    """

    sub_cigar = ''

    # finding the upstream flank in the read 
    upstream = ref.fetch(read.reference_name, locus_start - flank_length, locus_start)
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence, upstream)

    flank_position = 'upstream'
    compl_cigar = complement_cigar(sCigar)
    cindel_pos, cindel_type, cindel_length = get_indel_position(compl_cigar)
    cindel_pos += locus_start - flank_length + (query_begin - 1)
    if flank_position == relative_position:
        if indel_type == cindel_type and indel_position == cindel_pos and indel_length == cindel_length:
            return True

    # target_begin is the start of the alignment in the read
    # query_begin  is the start of the alignment in the reference flank
    # target_end   is the end of the alignment in the read
    # query_end    is the end of the alignment in the reference flank

    # if the alignment score is less than 90% of the flank length, return empty
    # making sure we identify the flanking sequence correctly in the read
    if alignment_score < int(2*(0.9*flank_length)): return [-1, -1, '']

    # reference position where the read (upstream flank) is starts to align
    read_reference_start = locus_start - flank_length + query_begin - 1

    # the read alignment to the flank; the position where the read aligns to the start of the flank
    if target_begin > 1: sub_cigar += f'{target_begin - 1}S'

    # because we are aligning the flank (ref sequence) to the read
    # here we complement the cigar to get the alignment of the read to the reference flank sequence
    sub_cigar += complement_cigar(sCigar)

    # if the flank aligns completely this should be equal to the repeat start
    uflank_reference_end = read_reference_start + query_end - query_begin + 1

    # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
    # these bases are at the end of the upstream flank
    uflank_remaining_length = flank_length - query_end
    uflank_read_end = target_end - 1
    # adding the remaining bases of the upstream flank as substitutions
    if uflank_remaining_length > 0:
        sub_cigar += f'{uflank_remaining_length}X'
        uflank_reference_end += uflank_remaining_length
        uflank_read_end += uflank_remaining_length

    # finding the downstream flank in the read
    downstream = ref.fetch(read.reference_name, locus_end, locus_end + flank_length)
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence, downstream)

    flank_position = 'downstream'
    compl_cigar = complement_cigar(sCigar)

    cindel_pos, cindel_type, cindel_length = get_indel_position(sCigar)
    cindel_pos += locus_end + (query_begin - 1)
    if flank_position == relative_position:
        if indel_type == cindel_type and indel_position == cindel_pos and indel_length == cindel_length:
            return True

    # target_begin is the start of the alignment in the read
    # query_begin  is the start of the alignment in the reference flank
    # target_end   is the end of the alignment in the read
    # query_end    is the end of the alignment in the reference flank

    # if the alignment score is less than 90% of the flank length, return empty
    # making sure we identify the flanking sequence correctly in the read
    if alignment_score < int(2*(0.9*flank_length)): return [-1, -1, '']

    dflank_cigar = ''

    # should be qual to the repeat end if the flank aligns completely
    dflank_reference_start = locus_end + query_begin - 1

    # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
    # these bases are at the start of downstream flank
    dflank_remaining_length = flank_length - query_end # number of bases in the flank that are not aligned
    dflank_read_start = target_begin - 1
    # adding the remaining bases of the downstream flank as substitutions
    if dflank_reference_start > locus_end:
        dflank_cigar += f'{dflank_remaining_length}X'     # adding the remaining bases as substitutions
        dflank_reference_start -= (dflank_remaining_length)
        dflank_read_start -= (dflank_remaining_length)

    dflank_cigar += complement_cigar(sCigar)

    # instead of adding the remaining downstream and upstream flank bases as substitutions
    # they could be included with the repeat sequence to be aligned with the reference repeat sequence
    # repeat_sequence = fasta.fetch(read.reference_name, uflank_read_end, dflank_read_start)
    # align the read repeat sequence with the reference repeat sequence

    # read_repeat_sequence = read.query_sequence[uflank_read_end:dflank_read_start]
    start_idx = uflank_read_end + 1
    end_idx = dflank_read_start
    allele_length = end_idx - start_idx
    if (start_idx >= end_idx): return [-1, -1, '']

    # aligning the sequence between the two flanks to the reference repeat
    alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(read.query_sequence[start_idx: end_idx],
                                                                                                   ref.fetch(read.reference_name, locus_start, locus_end))
    if  target_begin - query_begin > 0:  sub_cigar += f'{target_begin - query_begin}I'
    elif query_begin - target_begin > 0: sub_cigar += f'{query_begin - target_begin}D'

    if query_begin > 1: sub_cigar += f'{query_begin - 1}X'
    sub_cigar = complement_cigar(sCigar)
    remaining_subs = abs((locus_end-locus_start) - query_end)
    if remaining_subs > 0:
        sub_cigar += f'{remaining_subs}X'
        target_end += remaining_subs
        query_end += remaining_subs
    # if dflank_remaining_length > 0:
    #     sub_cigar += f'{dflank_remaining_length}D'
    if allele_length - target_end > 0:
        sub_cigar += f'{allele_length - target_end}I'

    # sub_cigar += dflank_cigar

    return [start_idx, end_idx, sub_cigar]


def update_homopolymer_coords(seq, locus_start, homopoly_positions):
    """
    Record all the homopolymer stretches of at least 3 bases within the repeat coordinates

    Args:
        seq (str): Reference sequence of the repeat region
        locus_start (int): Start position of the repeat region in the reference genome
        homopoly_positions (dict): Dictionary to store homopolymer positions
    
    Returns:
        None: Updates the homopoly_positions dictionary in place.
    """

    prev_N = seq[0]; start = -1
    for i, n in enumerate(seq[1:]):
        if n == prev_N:
            if start == -1: start = i
        else:
            if start != -1 and (i+1)-start >= 4:
                for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
                    homopoly_positions[c] = (i-start+1)-l
            start = -1
        prev_N = n

    if start != -1 and (i+1)-start >= 4:
        for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
            # for each position in the homopolymer stretch we record the length of the 
            # homopolymer nucleotides on the right
            homopoly_positions[c] = (i-start+1)-l


# def match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank, right_flank):
def match_jump(rpos, qpos, match_length, repeat_index, loci_keys, loci_coords, read_loci_info):
    """
    Records the matches and returns the last tracked repeat index
    Args:
        rpos (int): position in the reference genome where the match is occurring
        qpos (int): position in the read where the match is occurring
        match_length (int): length of the match
        repeat_index (int): index of the repeat locus in loci_coords
        loci_keys (list): list of keys for the loci
        loci_coords (list): list of coordinates for each locus
        read_loci_info (dict): dictionary containing information about each locus
    Returns:
        int: number of repeat indices to jump beyond the repeat where all positions are tracked
    """

    r = 0
    for r,coord in enumerate(loci_coords[repeat_index:]):
        ridx = repeat_index + r

        locus_key   = loci_keys[ridx]
        locus_start = loci_coords[ridx][0]
        locus_end   = loci_coords[ridx][1]

        # if the match has not crossed the repeat
        if rpos < locus_start: break
        # if the repeat is passed before the match
        if rpos - match_length > locus_end: continue

        # locus never tracked
        if not read_loci_info[locus_key]['tracked']:

            # updating the start position of the repeat on the read
            if locus_start <= rpos:
                read_loci_info[locus_key]['range'][0] = qpos - (rpos - locus_start)

            # updating the end position of the repeat on the read
            if locus_end <  rpos:
                read_loci_info[locus_key]['range'][1] = qpos - (rpos - locus_end)

            read_loci_info[locus_key]['tracked'] = True  # set tracked as true

        # if the locus is tracked but end position of the repeat in the read is not updated
        elif locus_end <= rpos:
            read_loci_info[locus_key]['range'][1] = qpos - (rpos - locus_end)

    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump


def deletion_jump(rpos, qpos, deletion_length, repeat_index, loci_keys, loci_coords,
                  read_loci_info, homopoly_positions):
    """
    Record the deletion in read and loci variations and return the last tracked repeat index
    Args:
        rpos (int): position in the reference genome where the deletion is occurring
        qpos (int): position in the read where the deletion is occurring
        deletion_length (int): length of the deletion
        repeat_index (int): index of the repeat locus in loci_coords
        loci_keys (list): list of keys for the loci
        loci_coords (list): list of coordinates for each locus
        read_loci_info (dict): dictionary to store variations in reads for each locus
        homopoly_positions (dict): dictionary containing positions of homopolymers
    Returns:
        int: number of repeat indices to jump beyond the repeat where all positions are tracked
    """

    # rpos - corresponds to the position in the reference after tracking the deletion
    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):
        ridx = repeat_index + r

        locus_key   = loci_keys[ridx]
        locus_start = coord[0]
        locus_end   = coord[1]

        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < locus_start: break

        # actual position in the reference where the deletion is occurring
        del_pos = rpos - deletion_length
        if del_pos > locus_end: continue

        if not read_loci_info[locus_key]['tracked']:
            # if the locus is not tracked and deletion is encountered beyond
            if locus_start <= rpos:
                read_loci_info[locus_key]['range'][0] = qpos
                read_loci_info[locus_key]['tracked'] = True    # set tracked as true

            # if the deletion crosses the end of the repeat
            if locus_end < rpos:
                read_loci_info[locus_key]['range'][1] = qpos

        # if the locus is tracked and the deletion crosses the end of the repeat
        elif locus_end < rpos:
            read_loci_info[locus_key]['range'][1] = qpos

        # getting the deletion length that covers only the repeat region
        locus_deletion_length = min(locus_end, rpos) - max(locus_start, del_pos)
        if (rpos >= locus_start) and (del_pos <= locus_end): # introduced to include length only if it comes inside repeat region
            if del_pos not in homopoly_positions:
                read_loci_info[locus_key]['ALEN']  -= locus_deletion_length
                read_loci_info[locus_key]['HALEN'] -= locus_deletion_length
            else:
                if locus_deletion_length <= homopoly_positions[del_pos]:
                    # if the deletion is only limited to the homopolymer positions
                    read_loci_info[locus_key]['HALEN'] -= locus_deletion_length
                else:
                    read_loci_info[locus_key]['ALEN']  -= locus_deletion_length
                    read_loci_info[locus_key]['HALEN'] -= locus_deletion_length

    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for coord in loci_coords[repeat_index:]:
            if coord[1] < rpos: jump += 1
            else: break

    return jump


def insertion_jump(rpos, qpos, insertion_length, insert, repeat_index, loci_keys, loci_coords, 
                   read_loci_info, read, ref, homopoly_positions, flank_length, aligned_start,
                   aligned_end, aligned_cigartuples):
    """
    Record the insertion in read and loci variations and return the last tracked repeat index
    Args:
        rpos (int): position in the reference genome where the insertion is occurring
        qpos (int): position in the read where the insertion is occurring
        insertion_length (int): length of the insertion
        insert (str): sequence of the inserted nucleotides
        repeat_index (int): index of the repeat locus in loci_coords
        loci_keys (list): list of keys for the loci
        loci_coords (list): list of coordinates for each locus
        read_loci_info (dict): dictionary to store variations in reads for each locus
        homopoly_positions (dict): dictionary containing positions of homopolymers
        flank_length (int): length of the flank to consider for insertion effects
    Returns:
        int: number of repeat indices to jump beyond the repeat where all positions are tracked
    """

    # a much shorter flank length within which if we detect an insertion/deletion we need to check if this affects the repeat

    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):

        locus_key = loci_keys[r+repeat_index]
        locus_start = coord[0]
        locus_end   = coord[1]

        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < locus_start - flank_length: break

        # if the insertion is happening beyond, the repeat in unaffected
        if rpos > locus_end + flank_length: continue

        if not read_loci_info[locus_key]['tracked']:
            # if the insertion is within the flank of the repeat
            if locus_start - flank_length <= rpos < locus_start:    
                # if an indel is found in the upstream flank of the repeat we realign the flank sequence to identify the coordinates of the repeat
                # result = detect_flank(read, ref, flank_length, locus_start, locus_end, 'upstream', 'I', insertion_length, rpos,
                #                       aligned_start, aligned_end, aligned_cigartuples)
                # # if realignment results in inertion at the same position we continue
                # if result == True or result[0] == -1 or result[1] == -1: pass
                # else:
                #     start_idx, end_idx, _ = result
                #     alen = end_idx - start_idx
                #     if alen > ((locus_end - locus_start)*10): pass
                #     else:
                #         read_loci_info[locus_key]['range'][0] = start_idx
                #         read_loci_info[locus_key]['range'][1] = end_idx
                #         read_loci_info[locus_key]['tracked']  = True
                #         read_loci_info[locus_key]['HALEN'] = end_idx - start_idx
                pass
           
            elif locus_end < rpos <= locus_end + flank_length:
                # if an indel is found in the downstream flank of the repeat we realign the flank sequence to identify the coordinates of the repeat
                # result = detect_flank(read, ref, flank_length, locus_start, locus_end, 'downstream', 'I', insertion_length, rpos,
                #                       aligned_start, aligned_end, aligned_cigartuples)
                # # if realignment results in insertion at the same position
                # if result == True or result[0] == -1 or result[1] == -1: pass
                # else:
                #     start_idx, end_idx, _ = result
                #     alen = end_idx - start_idx
                #     if alen > ((locus_end - locus_start)*10): pass
                #     else:
                #         read_loci_info[locus_key]['range'][0] = start_idx
                #         read_loci_info[locus_key]['range'][1] = end_idx
                #         read_loci_info[locus_key]['tracked']  = True
                #         read_loci_info[locus_key]['HALEN'] = end_idx - start_idx
                pass

            # if the locus is not tracked
            # what is the possible case that we reached in the middle of the repeat and it has not been tracked?
            elif locus_start <= rpos:
                if read_loci_info[locus_key]['range'][0] == -1:
                    read_loci_info[locus_key]['range'][0] = qpos-insertion_length
                read_loci_info[locus_key]['tracked']  = True    # set tracked as true

            # the insertion is happening at the end of the repeat
            if locus_end == rpos:
                if read_loci_info[locus_key]['range'][1] == -1:
                    read_loci_info[locus_key]['range'][1] = qpos
                # here jump can be done
        
        # if the locus is tracked
        else:
            # the insertion is exactly at the end position of the repeat
            if locus_end == rpos:
                if read_loci_info[locus_key]['range'][1] == -1:
                    read_loci_info[locus_key]['range'][1] = qpos

            # if the insertion is within the flank of the repeat
            elif locus_start - flank_length <= rpos and rpos < locus_start:
                # if an indel is found in the upstream flank of the repeat we realign the flank sequence to identify the coordinates of the repeat
                # result = detect_flank(read, ref, flank_length, locus_start, locus_end, 'upstream', 'I', insertion_length, rpos,
                #                       aligned_start, aligned_end, aligned_cigartuples)
                # # if realignment results in inertion at the same position we continue
                # if result == True or result[0] == -1 or result[1] == -1: pass
                # else:
                #     start_idx, end_idx, _ = result
                #     alen = end_idx - start_idx
                #     if alen > ((locus_end - locus_start)*10): pass
                #     else:
                #         read_loci_info[locus_key]['range'][0] = start_idx
                #         read_loci_info[locus_key]['range'][1] = end_idx
                #         read_loci_info[locus_key]['tracked']  = True
                #         read_loci_info[locus_key]['HALEN'] = end_idx - start_idx
                pass

            elif locus_end < rpos and rpos <= locus_end + flank_length:
                # if an indel is found in the downstream flank of the repeat we realign the flank sequence to identify the coordinates of the repeat
                # result = detect_flank(read, ref, flank_length, locus_start, locus_end, 'downstream', 'I', insertion_length, rpos,
                #                       aligned_start, aligned_end, aligned_cigartuples)
                # # if realignment results in insertion at the same position
                # if result == True or result[0] == -1 or result[1] == -1: pass
                # else:
                #     start_idx, end_idx, _ = result
                #     alen = end_idx - start_idx
                #     if alen > ((locus_end - locus_start)*10): pass
                #     else:
                #         read_loci_info[locus_key]['range'][0] = start_idx
                #         read_loci_info[locus_key]['range'][1] = end_idx
                #         read_loci_info[locus_key]['tracked']  = True
                #         read_loci_info[locus_key]['HALEN'] = end_idx - start_idx
                pass

        # if the insert falls in the flank of the repeat we then try to find the flanks by sequence
        if locus_start <= rpos <= locus_end: # introduced to include length only if it comes inside repeat region
            if len(set(insert)) == 1:
                # only if the insertion is a homopolymer; consider it as homopolymer insertion
                if read_loci_info[locus_key]['range'][0] == -1 or read_loci_info[locus_key]['range'][1] == -1:
                    read_loci_info[locus_key]['HALEN'] += insertion_length
            else:
                if read_loci_info[locus_key]['range'][0] == -1 or read_loci_info[locus_key]['range'][1] == -1:
                    read_loci_info[locus_key]['ALEN']  += insertion_length
                    read_loci_info[locus_key]['HALEN'] += insertion_length

    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump
