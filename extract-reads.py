import pysam
import argparse
from Bio.Align import PairwiseAligner

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


def parse_cigar(cigar_tuples, read_start, repeat_start, repeat_end, ref_sequence, query_sequence):
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
        if cigar[0] == 4:
           # soft clipped - these bases are part of the read sequence but do not
           #                affect the reference position.
           qpos += cigar[1] 

        elif cigar[0] == 2:     # deletion
            deletion_length = cigar[1]

            for i in range(deletion_length):
                # iterating over the each deleted base
                if rpos + i == repeat_start and start_idx == False:
                    # if the deletion covers the repeat start
                    start_idx = qpos

                if rpos + i == repeat_end:
                    # deletion covers the repeat end
                    end_idx = qpos

                if start_idx != False and end_idx == False:
                    # if the position is within the repeat
                    sub_cigar += 'D'
            # move the reference position
            rpos += deletion_length

        elif cigar[0] == 1:     # insertion
            insert_length = cigar[1]

            if rpos == repeat_start and start_idx == False:
                # if the insert is before the repeat include the sequence within the repeat
                start_idx = qpos
            for i in range(insert_length):
                if start_idx != False and end_idx == False:
                    # if insert within the repeat
                    sub_cigar += 'I'
            # move the query position
            qpos += insert_length

        elif cigar[0] == 0: # match (both equals & difference)
            match_len = cigar[1]

            for i in range(match_len):
                # iterating over the match positions
                if rpos + i == repeat_start and start_idx == False:
                    # if match covers the repeat start
                    start_idx = qpos + i

                if rpos + i == repeat_end:
                    # if match covers the repeat end
                    end_idx = qpos + i

                if start_idx != False and end_idx == False:
                    # if match within the repeat
                    if query_sequence[qpos+i] == ref_sequence[rpos+i-read_start]:
                        # if the bases match
                        sub_cigar += 'M'
                    else:
                        # if its a mismatch
                        sub_cigar += 'X'
            # move both reference and query positions
            rpos += match_len; qpos += match_len

        if rpos > repeat_end:
            if end_idx == False:
                # if the end of the repeat is not covered by the read
                # but the read has moved beyond the repeat
                end_idx = qpos - (rpos - repeat_end)
            # if position moved beyond repeat
            return [start_idx, end_idx, sub_cigar]

    if end_idx == False:
                # if the end of the repeat is not covered by the read
                # but the read has moved beyond the repeat
                end_idx = qpos - (rpos - repeat_end)
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
            motif_len = int(line[4])

            if chrom != 'chrX': continue

            reference_repseq = ""
            for bam in bams:
                read_data = []
                if chrom not in bam.references: continue
                reads = bam.fetch(chrom, repeat_start, repeat_end)
                check = False
                for read in reads:
                    if read.reference_start < repeat_start-10 and read.reference_end > repeat_end + 10:
                        print("Aligned in repeat:", read.reference_start, read.reference_end, read.cigarstring, sep='\t')
                        print(read.query_sequence, sep='\t')
                        if not read.has_tag('MD'): continue
                        if not check:
                            reference_repseq = read.get_reference_sequence()[(repeat_start - read.reference_start): (repeat_start - read.reference_start) + repeat_end-repeat_start]
                            # print(chrom, repeat_start, repeat_end, reference_repseq, cigar, sep='\t')
                            print(*line, sep='\t')
                            check = True
                        
                        print(read.has_tag('MD'))
                        start_idx, end_idx, sub_cigar = parse_cigar(read.cigartuples, read.reference_start, repeat_start, repeat_end,
                                                                    read.get_reference_sequence(), read.query_sequence)
                        
                        read_repseq = read.query_sequence[start_idx: end_idx]
                        motif = reference_repseq[:motif_len]
                        # RD2RP_CIGAR, RF2RD2RP_CIGAR, tags = convert_CIGAR(cigar, sub_cigar, motif, motif_len, read_repseq)

                        read_data.append([chrom, repeat_start, repeat_end, read.query_name, start_idx, end_idx])

                    else:
                        if read.cigartuples[0][0] == 4:
                            # soft clipped - these bases are part of the read sequence but do not
                            #                affect the reference position.
                            soft_clipped = read.query_sequence[:read.cigartuples[0][1]]
                            sclip_len = read.cigartuples[0][1]
                            if read.has_tag('MD'):
                                # if the read has a MD tag
                                upstream = fasta.fetch(chrom, repeat_start - 50, repeat_start)
                                aligner = PairwiseAligner()
                                aligner.mode = 'local'
                                best_alingment_score = 0
                                best_position = 0
                                best_alignment = None
                                for i in range(len(soft_clipped)):
                                    alignment = aligner.align(soft_clipped[i:i+50], upstream)
                                    if alignment[0].score > best_alingment_score and alignment[0].score > 40:
                                        best_alingment_score = alignment[0].score
                                        best_position = i
                                        best_alignment = alignment[0]
                                start_idx, end_idx, sub_cigar = parse_cigar([(4, best_position), (0, 50), (1,sclip_len-best_position-50)] + read.cigartuples[1:], read.reference_start, repeat_start, repeat_end,
                                                                            read.get_reference_sequence(), read.query_sequence)
                                read_data.append([chrom, repeat_start, repeat_end, read.query_name, start_idx, end_idx])
                            else: print("No MD tag")

                read_data = sorted(read_data, key=lambda x: x[2])
                for data in read_data:
                    # qualities = ':'.join([str(x) for x in data[-1]])
                    # print(*data[:-1], qualities, sep='\t')
                    print(*data, sep='\t')

    for bam in bams: bam.close()
    fasta.close()


if __name__ == "__main__":
    args = parse_args()

    aln_format = 'rb'
    if args.aln_format == 'cram': aln_format = 'rc'
    elif args.aln_format == 'sam': aln_format = 'r'

    extract_reads(args.bed, args.bam, args.ref_fasta, aln_format)
