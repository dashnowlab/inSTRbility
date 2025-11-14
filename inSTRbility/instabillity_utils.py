from statistics import mode
from lib_ssw.pyssw import align_pair
from collections import Counter
import numpy as np
import sys

import pyabpoa as pa


def consensus_seq_poa(seqs):
    """
    Generate a consensus sequence using the Partial Order Alignment (POA) algorithm.
    Args:
        seqs (list): List of sequences to be aligned.
    Returns:
        str: The consensus sequence generated from the input sequences.
    """

    if len(seqs) < 7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    abpoa  = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]


def nonatomic(motif, motif_length):
    """
    Check if the sequence is composed of pure stretches of a given motif length.
    Args:
        sequence (str): The input sequence to check.
        motif_length (int): The length of the motif to check for.
    Returns:
        bool: True if the sequence is composed of pure stretches, False otherwise.
    """

    for i in range(1, motif_length//2 + 1):
        if motif_length % i != 0: continue
        amotif = motif[:i]; units = motif_length // i
        if amotif * units == motif:
            return True
    return False


def get_pure_stretches(sequence, mlen):
    """
    Find the pure stretches of a given motif length in a sequence.
    Args:
        sequence (str): The input sequence to search for repeats.
        mlen (int): The length of the motif to search for.
    Returns:
        list: A list of tuples containing the motif, start index, and end index of each repeat.
    """

    repeats = []
    slen = len(sequence)
    i = 0
    while i < slen - mlen + 1:
        motif = sequence[i:i + mlen]

        if motif == sequence[i+mlen: i+2*mlen]:
            repeat_start = i
            repeat_end = i + 2*mlen
            m = 0
            j = i + 2 * mlen
            while j < slen and sequence[j] == motif[m]:
                m += 1
                j += 1
                repeat_end = j
                if m == mlen: m = 0
            if not nonatomic(motif, mlen):
                repeats.append((motif, repeat_start, repeat_end))

            i = j - (mlen - 1)

        i += 1

    return repeats


def extend_cigar(target_begin, target_end, query_begin, query_end, cigar, consensus_seq, read_seq):
    """
    Adjust the CIGAR string to account for the complete alignment of the target and query sequence including the soft-clipped regions.

    Args:
        target_begin (int): Start position of the target sequence.
        target_end (int): End position of the target sequence.
        query_begin (int): Start position of the query sequence.
        query_end (int): End position of the query sequence.
        cigar (str): The original CIGAR string.
        consensus_seq (str): The consensus sequence.
        read_seq (str): The read sequence.
    Returns:
        str: The adjusted CIGAR string.
    """

    target_length = len(consensus_seq)
    query_length = len(read_seq)

    target_begin -= 1; query_begin -= 1

    if target_begin > 0:
        if query_begin == 0:
            cigar = f'{target_begin}D' + cigar
        elif query_begin > 0:
            if query_begin < target_begin:
                cigar = f'{target_begin-query_begin}D' + f'{query_begin}X' + cigar
            elif query_begin > target_begin:
                cigar = f'{query_begin-target_begin}I' + f'{target_begin}X' + cigar

    elif target_begin == 0:
        if query_begin > 0:
            cigar = f'{query_begin}I' + cigar

    if target_end < target_length:
        if query_end == query_length:
            cigar += f'{target_length - target_end}D'
        elif query_end < query_length:
            target_missed = target_length - target_end
            query_missed = query_length - query_end
            if target_missed > query_missed:
                cigar += f'{query_missed}X' + f'{target_missed-query_missed}D'
            elif target_missed == query_missed:
                cigar += f'{target_missed}X'
            elif target_missed < query_missed:
                cigar += f'{target_missed}X' + f'{query_missed-target_missed}I'

    elif target_end == target_length:
        if query_end < query_length:
            cigar += f'{query_length - query_end}I'

    return cigar



def check_insert(insert, motif):
    """
    Check if the insertion is a pure stretch of the motif.

    Args:
        insert (str): The insertion sequence.
        motif (str): The motif to check against.

    Returns:
        bool: True if the insertion is a pure stretch of the motif, False otherwise.
    """

    ref = motif * ((len(insert) // len(motif) + 2))
    score, orientation, target_begin, target_end, query_begin, query_end, cigar = align_pair(ref, insert)
    return score >= (0.9 * 2 * len(insert))


def pure_stretch_variations(cigar, query, pure_stretches):
    """
    Process the CIGAR string to calculate the variations in lengths of pure stretches.
    Args:
        cigar (str): The CIGAR string representing the alignment.
        pure_stretches (list): List of tuples containing motif, start, and end positions of pure stretches.
    Returns:
        list: A list of variation in lengths of the pure stretches.
    """

    qpos = 0; rpos = 0
    clen = ''; ctype = ''
    alens = [0 for a in pure_stretches]
    for c in cigar:
        if c.isdigit(): clen += c

        else:
            ctype = c
            clen = int(clen) if clen else 0
            if ctype == 'D':
                deleted_positions = set(range(rpos, rpos + clen))
                for r, repeat in enumerate(pure_stretches):
                    motif, start, end = repeat
                    repeat_positions = set(range(start, end))
                    if deleted_positions & repeat_positions:
                        alens[r] -= len(deleted_positions & repeat_positions)
                rpos += clen

            elif ctype == 'I':
                insert = query[qpos:qpos + clen]
                for r, repeat in enumerate(pure_stretches):
                    motif, start, end = repeat

                    insert_check = False
                    if motif == 'random': insert_check = True
                    else: insert_check = check_insert(insert, motif)
                    if qpos >= start and qpos <= end and insert_check:
                        alens[r] += clen
                qpos += clen

            else:
                rpos += clen; qpos += clen
            clen = ''

    return alens


def consensus_allele(locus_key, LOCI_VARS, LOCI_INFO, motif_length, pure_stretch, num_hpreads, out):
    """
    Generate a consensus allele for each locus based on the read indices.

    Args:
        locus_key (str): The key representing the locus.
        LOCI_VARS (dict): Dictionary containing locus variables.
        LOCI_INFO (dict): Dictionary containing locus information.
        motif_length (int): The length of the motif to consider for pure stretches.
        pure_stretch (bool): Whether to consider pure stretches in the analysis.

    Returns:
        dict: A dictionary with locus keys and their corresponding consensus alleles.
    """

    locus_chrom = LOCI_INFO[locus_key][0]
    locus_start = int(LOCI_INFO[locus_key][1])
    locus_end   = int(LOCI_INFO[locus_key][2])
    read_indices = LOCI_VARS[locus_key]['RDS']
    sorted_ridx  = sorted(list(range(0, len(read_indices))), key=lambda x: LOCI_VARS[locus_key]['HP'][x] if LOCI_VARS[locus_key]['HP'][x] is not None else 3)

    haplogrouped_reads = {}
    haplogrouped_sequences = {}
    for ridx in sorted_ridx:
        read_idx = LOCI_VARS[locus_key]['RDS'][ridx]
        fail_messages = {1: 'low-reads', 2: 'nosig-snps', 3: 'low-read-snps', 4: 'low-reads-cluster', 6: 'low-reads', 10: 'pass'}
        haplogroup = LOCI_VARS[locus_key]['HP'][ridx]

        if haplogroup is None: haplogroup = 0
        if haplogroup not in haplogrouped_reads:
            haplogrouped_reads[haplogroup] = []
            haplogrouped_sequences[haplogroup] = []

        haplogrouped_reads[haplogroup].append(read_idx)
        haplogrouped_sequences[haplogroup].append(LOCI_VARS[locus_key]['S'][read_idx])

    # filtering out the haplogroups with low read support
    removed_haplogroups = []
    for haplogroup, sequences in haplogrouped_sequences.items():
        haplogroup_alen = mode([len(seq) for seq in sequences])

        # if either the haplogroup has less number of reads supporting it or length of the haplogroup is less than 12, remove it
        if len(sequences) >= num_hpreads or haplogroup_alen >= 12: continue
        removed_haplogroups.append(haplogroup)

    for haplogroup in removed_haplogroups:
        del haplogrouped_reads[haplogroup]
        del haplogrouped_sequences[haplogroup]

    # if there are less than 2 haplogroups, we cannot proceed
    if 1 not in haplogrouped_reads or 2 not in haplogrouped_reads: return

    purestretch_data = {}

    for haplogroup, sequences in haplogrouped_sequences.items():
        haplogroup_alen = mode([len(seq) for seq in sequences])
        if haplogroup_alen < 12:
            if int(haplogroup) > 0: return
            continue

        consensus_seq = consensus_seq_poa([seq for seq in sequences if len(seq) == haplogroup_alen])

        # for non-repetitive loci, we do not identify pure stretches
        if motif_length > 0:
            if pure_stretch:
                pure_stretches = get_pure_stretches(consensus_seq, motif_length)
            else:
                pure_stretches = [(LOCI_INFO[locus_key][3], 0, len(consensus_seq))]

        else: pure_stretches = [('random', 0, len(consensus_seq))]

        if not pure_stretches:
            if int(haplogroup) > 0: return
            continue

        # stores the length variation for each pure stretch in each read
        read_lenvars = [[] for _ in pure_stretches]

        for read_idx in haplogrouped_reads[haplogroup]:
            read_seq = LOCI_VARS[locus_key]['S'][read_idx]
            if len(read_seq) == 0: continue
            if len(consensus_seq) == len(read_seq):
                for _ in range(len(pure_stretches)):
                    read_lenvars[_].append(0)
                continue

            score, orientation, target_begin, target_end, query_begin, query_end, cigar = align_pair(consensus_seq, read_seq)
            new_cigar = extend_cigar(target_begin, target_end, query_begin, query_end, cigar, consensus_seq, read_seq)
            purestretch_lenvars = pure_stretch_variations(new_cigar, read_seq, pure_stretches)
            for i, var in enumerate(purestretch_lenvars):
                read_lenvars[i].append(var)

        for i in range(len(pure_stretches)):
            locus = pure_stretches[i]
            plocus_key = f'{locus_chrom}:{locus_start+locus[1]}-{locus_start+locus[2]}'

            motif = locus[0]
            alen = locus[2] - locus[1]
            mad = round(np.mean([abs(x) for x in read_lenvars[i]]), 3)
            filtered_vars = [x for x in read_lenvars[i] if x > 0]
            if not filtered_vars:
                mad2 = 0.0
            else:
                mad2 = round(np.median([abs(x) for x in filtered_vars]), 3)

            # if the length of the allele is less that 3 times the motif length, skip it
            if alen < 3*motif_length or alen < 12: continue

            allele_counts = Counter(read_lenvars[i])
            allele_counts = ';'.join([f'{x[0]}:{x[1]}' for x in allele_counts.most_common()])

            if plocus_key not in purestretch_data:
                purestretch_data[plocus_key] = {haplogroup: {}}
            elif haplogroup not in purestretch_data[plocus_key]:
                purestretch_data[plocus_key][haplogroup] = {}

            purestretch_data[plocus_key][haplogroup]['alen'] = alen
            purestretch_data[plocus_key][haplogroup]['motif'] = motif
            purestretch_data[plocus_key][haplogroup]['mad'] = mad
            purestretch_data[plocus_key][haplogroup]['mad2'] = mad2
            purestretch_data[plocus_key][haplogroup]['allele_counts'] = allele_counts

    for plocus_key, groups in purestretch_data.items():
        if 1 not in groups or 2 not in groups: continue
        for haplogroup, data in purestretch_data[plocus_key].items():
            motif = data['motif']
            alen = data['alen']
            mad = data['mad']
            mad2 = data['mad2']
            allele_counts = data['allele_counts']

            print(plocus_key, haplogroup, motif, alen, mad, mad2, allele_counts, sep='\t', file=out)
