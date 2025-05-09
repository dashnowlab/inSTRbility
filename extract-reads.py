import pysam
import argparse
from lib_ssw.pyssw import align_pair
import statistics as stats

"""
NOTE: Though samtools uses a position with a 1 based coordinate system. read.reference_start returns the position w.r.t 
      a 0-based coordinate system.
"""

def parse_args():
    parser = argparse.ArgumentParser(prog='extract-reads', description="Extract the reads aligned at a locus in a BAM file and analyse.")

    parser.add_argument('-bam', required=True, nargs="+", type=str, dest="bam", help="Input BAM file from which the reads to be extracted")
    parser.add_argument('-bed', required=True, type=str, dest="bed", help="Input regions file for which the reads are extracted")
    parser.add_argument('-ref',  required=True, type=str, dest="ref_fasta", help="The reference fasta genome file")

    parser.add_argument('--aln-format', default='bam', type=str, dest="aln_format", help="Format of the alignment files. Choose from bam/cram/sam. Default: bam")

    args = parser.parse_args()

    return args


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


def parse_cigar(cigar_tuples, read_start, repeat_start, repeat_end):
    """
    Parses read alignment CIGAR and returns coordinates of the repeat within the read and the CIGAR
    corresponding to the repeat sequence.

    Args:
        cigar_tuples:   cigar tuples from pysam record
        read_start:     reference position of start of the read alignment
        repeat_start:   start coordinate of repeat in reference
        repeat_end:     end coordinate of repeat in reference
        ref_sequence:   reference sequence to which read is aligned
        query_sequence: sequence of the read

    Returns:
        start and end coordinates of read sequence aligning to the repeat region
        CIGAR string of alignment within repeat sequence
        [start, end, sub_cigar]
    """
    rpos = read_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0            # starts from 0 the sub string the read sequence in python

    # the read start and the repeat start are both on a 0 based coordinate system

    start_idx = False; end_idx = False
    sub_cigar = ''

    for c, cigar in enumerate(cigar_tuples):
        # print(cigar, rpos, qpos, start_idx, end_idx, repeat_start, repeat_end, sub_cigar, sep='\t')
        if cigar[0] == 4:
           # soft clipped - these bases are part of the read sequence but do not
           #                affect the reference position.
           qpos += cigar[1] 

        elif cigar[0] == 2:     # deletion
            deletion_length = cigar[1]
            if start_idx == False and rpos + deletion_length >= repeat_start:
                start_idx = qpos
                if rpos + deletion_length >= repeat_end:
                    end_idx = qpos
                    sub_cigar += f'{repeat_end - repeat_start}D'
                else: sub_cigar += f'{rpos + deletion_length - repeat_start}D'
            
            elif start_idx != False and end_idx == False:
                if rpos + deletion_length >= repeat_end:
                    end_idx = qpos
                    sub_cigar += f'{repeat_end-rpos}D'
                else: sub_cigar += f'{deletion_length}D'

            # move the reference position
            rpos += deletion_length

        elif cigar[0] == 1:     # insertion
            insert_length = cigar[1]

            if rpos == repeat_start and start_idx == False:
                # if the insert is before the repeat include the sequence within the repeat
                start_idx = qpos
                sub_cigar += f'{insert_length}I'
            elif start_idx != False and end_idx == False:
                sub_cigar += f'{insert_length}I'
            elif start_idx != False and rpos == repeat_end:
                sub_cigar += f'{insert_length}I'
                end_idx = qpos + insert_length
            # move the query position
            qpos += insert_length

        elif cigar[0] == 0 or cigar[0] == 7 or cigar[0] == 8: # match (both equals & difference)
            match_len = cigar[1]
            ctype = 'X' if cigar[0] == 8 else 'M'

            if start_idx == False and rpos + match_len > repeat_start:
                start_idx = qpos + (repeat_start - rpos)
                if rpos + match_len >= repeat_end:
                    end_idx = qpos + (repeat_end - rpos) 
                    sub_cigar += f'{repeat_end - repeat_start}{ctype}'
                else: sub_cigar += f'{rpos + match_len - repeat_start}{ctype}'
            
            elif start_idx != False and end_idx == False:
                if rpos + match_len >= repeat_end:
                    end_idx = qpos + (repeat_end - rpos)
                    sub_cigar += f'{repeat_end - rpos}{ctype}'
                else:
                    sub_cigar += f'{match_len}{ctype}'

            # move both reference and query positions
            rpos += match_len; qpos += match_len

        if rpos > repeat_end:
            # if position moved beyond repeat
            return [start_idx, end_idx, sub_cigar]

    return [start_idx, end_idx, sub_cigar]


def extract_reads(bed_file, bam_files, ref_fasta, aln_format):
    """
    Given the set of repeat loci. The function extracts reads aligning to the repeat locus
    within a given set of BAMs

    Args:
        bed_file:   list of input of repeat regions as a BED file
        bam_files:  list of input BAM files (Ribbit input)
        aln_format: alignment file format (bam/sam/cram)
    """
    bams = [pysam.AlignmentFile(bam_file, aln_format, check_sq=False) for bam_file in bam_files]
    fasta = pysam.FastaFile(ref_fasta)

    with open(bed_file) as fh:
        for line in fh:
            line = line.strip().split('\t')
            chrom = line[0]
            repeat_start = int(line[1]); repeat_end = int(line[2])
            cigar = line[-1]
            motif_len = len(line[3])
            start_idx = -1; end_idx = -1

            if chrom != 'chrX': continue

            for bam in bams:
                # Get the file name
                bam_id = bam.filename.decode('utf-8').split('.')[0]
                read_data = []
                if chrom not in bam.references: continue
                reads = bam.fetch(chrom, repeat_start, repeat_end)
                check = False
                for read in reads:

                    if read.reference_start < repeat_start-10 and read.reference_end > repeat_end + 10:
                        start_idx, end_idx, sub_cigar = parse_cigar(read.cigartuples, read.reference_start, repeat_start, repeat_end)

                    else:
                        flank_len = 50
                        startclip_len = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
                        endclip_len = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
                        include_startclip = False; include_endclip = False
                        sclip_cigartuples = []; eclip_cigartuples = []; new_cigartuples = []
                        new_reference_start = read.reference_start

                        if read.reference_start >= repeat_start - 10 and read.reference_start - startclip_len < repeat_start - 10:
                            include_startclip = True
                        if read.reference_end <= repeat_end + 10 and read.reference_end + endclip_len > repeat_end + 10:
                            include_endclip = True
                        if include_startclip:
                            soft_clipped = read.query_sequence[:startclip_len]
                            if read.reference_start - startclip_len >= repeat_start - 10: continue

                            upstream = fasta.fetch(chrom, read.reference_start - flank_len, read.reference_start)
                            alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(soft_clipped, upstream)
                            
                            sclip_cigartuples.append((4, target_begin))
                            flank_reference_start = (read.reference_start - flank_len) + (query_begin - 1)
                            sclip_cigartuples += convert_cigar(sCigar)
                            flank_reference_end = flank_reference_start + query_end
                            sclip_insert_len = startclip_len - target_end
                            
                            if read.reference_start - flank_reference_end > 0:
                                sclip_cigartuples.append((0, read.reference_start - flank_reference_end))
                            sclip_cigartuples.append((1, sclip_insert_len - (read.reference_start - flank_reference_end)))
                            new_reference_start = flank_reference_start

                        
                        if include_endclip:
                            soft_clipped = read.query_sequence[-endclip_len:]
                            if read.reference_end + endclip_len <= repeat_end + 10: continue

                            downstream = fasta.fetch(chrom, read.reference_end, read.reference_end + flank_len)
                            alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar = align_pair(soft_clipped, downstream)

                            eclip_insert_len = target_begin - 1 - query_begin - 1
                            if eclip_insert_len > 0:
                                eclip_cigartuples.append((1, eclip_insert_len))
                            if query_begin - 1 > 0: eclip_cigartuples.append((8, query_begin - 1))
                            eclip_cigartuples += convert_cigar(sCigar)
                            if endclip_len - target_end > 0: eclip_cigartuples.append((4, endclip_len - target_end))
                        
                        if len(sclip_cigartuples) > 0:
                            new_cigartuples = sclip_cigartuples + read.cigartuples[1:]
                        else: new_cigartuples = read.cigartuples
                        if len(eclip_cigartuples) > 0:
                            new_cigartuples = new_cigartuples[:-1] + eclip_cigartuples
                        
                        # print(read.cigarstring, convert_cigartuples(new_cigartuples))
                        start_idx, end_idx, sub_cigar = parse_cigar(new_cigartuples, new_reference_start, repeat_start, repeat_end)


                    # Get methylation
                    chunk_meth = []
                    mods = read.modified_bases_forward
                    for modtype in mods:
                        if modtype[0] == 'C' and modtype[2] == 'm':
                            for pos, qual in mods[modtype]:
                                if pos >= start_idx and pos < end_idx: # Need to check the position logic here
                                    prob = qual/256
                                    if qual > 0:
                                        chunk_meth.append(prob)
                    # median meth
                    if len(chunk_meth) > 0:
                        med_meth = round(stats.median(chunk_meth), 3)
                    else:
                        med_meth = None
                    allele_len = round((end_idx - start_idx)/motif_len, 2)

                    # print(read.cigarstring, read.reference_start)
                    read_data.append([f"{chrom}:{repeat_start}-{repeat_end}", read.query_name, start_idx, end_idx, sub_cigar, allele_len, med_meth])
                    print(*read_data[-1], sep='\t')

                # read_data = sorted(read_data, key=lambda x: x[5])
                # Print header
                # print(f"#sample\tlocus\tread_name\tread_repeat_start\tread_repeat_end\trepeat_cigar\tallele_length\tmedian_meth")
                # for data in read_data:
                #     print(*data, sep='\t')

    for bam in bams: bam.close()
    fasta.close()


if __name__ == "__main__":
    args = parse_args()

    aln_format = 'rb'
    if args.aln_format == 'cram': aln_format = 'rc'
    elif args.aln_format == 'sam': aln_format = 'r'

    extract_reads(args.bed, args.bam, args.ref_fasta, aln_format)
