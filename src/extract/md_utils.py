import bisect

def update_global_snps(ref_start, pos, READ_VARS, SNPS, read_index, read_sequence, read_basequalities,
                       SORT_SNPS, insertion_point, qpos):

    rpos = ref_start+pos
    for ins in insertion_point:
        if ins<rpos:
            qpos+=insertion_point[ins]
        elif ins>rpos: break

    base_phredQ = read_basequalities[qpos]
    sub_char = read_sequence[qpos]
    READ_VARS[read_index]['X'].add(rpos)
    if rpos not in SNPS:
        SNPS[rpos] = { 'COV': 1, sub_char: {read_index}, 'Q': {read_index : base_phredQ} }
        bisect.insort(SORT_SNPS, rpos)
    else:
        SNPS[rpos]['COV'] += 1
        SNPS[rpos]['Q'][read_index] = base_phredQ
        if sub_char in SNPS[rpos]: 
            SNPS[rpos][sub_char].add(read_index)

        else: SNPS[rpos][sub_char] = {read_index}


def parse_mdtag(MD_tag, qpos, ref_start, READ_VARS, SNPS, read_index, read_basequalities,
                read_sequence, SORT_SNPS, insertion_point):

    if SORT_SNPS == None: SORT_SNPS = []

    base = 0
    sub_base='0'
    sub_char=''
    pos=0

    deletion = False
    replacing = False

    for i in MD_tag:

        if deletion:
            if i.isalpha():
                base += 1
                continue
            else: deletion = False


        if i.isnumeric():
            sub_base += i
            if sub_char != '':
                update_global_snps(ref_start, pos, READ_VARS, SNPS, read_index, read_sequence, read_basequalities,
                                   SORT_SNPS, insertion_point, qpos)
                replacing = False
                qpos += 1
                sub_char = ''

        elif i.isalpha():
            replacing = True

            if sub_char == '':
                base += int(sub_base) + 1
                pos = base - 1
                qpos += int(sub_base)
            else:
                base += 1
                update_global_snps(ref_start, pos, READ_VARS, SNPS, read_index, read_sequence, read_basequalities,
                                     SORT_SNPS, insertion_point, qpos)
                pos = base - 1
                qpos += 1
                sub_char = ''

            sub_base = ''
            sub_char += i

        else: #i == '^':
            deletion = True

            if replacing:
                if sub_char != '':
                    update_global_snps(ref_start, pos, READ_VARS, SNPS, read_index, read_sequence, read_basequalities,
                                         SORT_SNPS, insertion_point, qpos)
                    replacing = False
                    qpos += 1
                    sub_char = ''
            else:
                base += int(sub_base)
                qpos += int(sub_base)
                sub_base = ''

    if sub_char != '':
        update_global_snps(ref_start, pos, READ_VARS, SNPS, read_index, read_sequence, read_basequalities,
                             SORT_SNPS, insertion_point, qpos)