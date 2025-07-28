#!/usr/bin/env python

"""
    inSTRbility is a tool designed to analyse somatic instability at tandem repeat loci
"""

import sys, os
import pysam
import timeit
import argparse as ap
from multiprocessing import Process

from version import __version__
from process_reads import *

def parse_args():
    """
    Parse command line arguments.
    """
    parser = ap.ArgumentParser(prog='inSTRbility', description="inSTRbility is a program that let\'s you analyse the instability at tandem repeat loci.")
    parser._action_groups.pop()

    print("inSTRbility - Analysing somatic instability at tandem repeat loci.\nDashnow Lab\n")

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-ref', required=True, type=str, dest='ref', metavar='<FILE>', help="The reference fasta genome file.")
    required.add_argument('-bam', required=True, nargs='+', type=str, dest='bam', metavar='<FILE>', help="Input alignment files.")
    required.add_argument('-bed', '--regions', dest='bed', required=True, metavar='<FILE>', help="Input regions file in bed format.")

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--aln-format', type=str, metavar='<STR>', default='bam', help="Format of the alignment files. Choose from bam/cram/sam. Default: bam")

    optional.add_argument('-o',  '--output', dest='output', type=str, metavar='<FILE>', default='', help='Name of the output file.')

    optional.add_argument('--contigs', nargs='+', help='Contigs to get genotyped [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be genotyped.')
    optional.add_argument('--flank',   type=int,  metavar='<INT>', default=30,  help='Length of the flanking region (in base pairs) to search for indels with a repeat in it. [default: 10]')
    optional.add_argument('-q', '--base-qual', type=int,  metavar='<INT>', dest='base_qual', default=10,   help='Minimum average base call quality of the whole read and also the bases corresponding to the repeat region. [default: 10]')
    optional.add_argument('--map-qual', type=int,  metavar='<INT>', dest='mapq', default=5,   help='Minimum mapping quality of the reads to be considered. [default: 5]')
    optional.add_argument('--min-reads',      type=int,  metavar='<INT>', default=10,  help='Minimum read coverage after quality cutoff at a locus to be genotyped. [default: 10]')
    optional.add_argument('--max-reads',      type=int,  metavar='<INT>', default=300, help='Maximum number of reads to be used for genotyping a locus. [default: 100]')

    optional.add_argument('--snp-dist',  type=int,   metavar='<INT>',   default=10000, help='Maximum distance of the SNP from repeat region to be considered for phasing. [default: 10000]')
    optional.add_argument('--snp-count', type=int,   metavar='<INT>',   default=3,    help='Number of SNPs to be considered for phasing (minimum value = 1). [default: 3]')
    optional.add_argument('--snp-qual',  type=int,   metavar='<INT>',   default=13,   help='Minimum basecall quality at the SNP position to be considered for phasing. [default: 13]')
    optional.add_argument('--snp-reads', type=float, metavar='<FLOAT>', default=0.2,  help='Minimum fraction of reads supporting an SNP to be used for phasing. [default: 0.25]')

    optional.add_argument('--fraction-phased-reads', type=float, metavar='<FLOAT>', default=0.4, help='Minimum fraction of total read contribution from the phased read clusters. [default: 0.4]')

    optional.add_argument('-t',  '--threads', type=int, metavar='<INT>', default=1, help='number of processor. [default: 1]')
    optional.add_argument('--haplotag', type=str, metavar='<STR>', default=None, help='use haplotagged information for phasing. eg: [HP]. [default: None]')
    optional.add_argument('--karyotype', nargs='+', help='karyotype of the samples [XY XX]')

    optional.add_argument('-v',   '--version', action='version',    version=f'inSTRbility version {__version__}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


def fasta_check(path):
    """
    Check if the provided FASTA file is valid.

    Args:
        path (str): Path to the FASTA file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid FASTA file.
        OSError: If there is an error reading the file.
    """
    try:
        f = pysam.FastaFile(path)
        f.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid FASTA file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()


def bam_check(path, aln_format):
    """
    Check if the provided BAM file is valid and sorted by coordinate.

    Args:
        path (str): Path to the BAM file.
        aln_format (str): Format of the alignment file.

    Raises:
        FileNotFoundError: If the file does not exist.
        SortOrderError: If the BAM file is not sorted by coordinate.
        ValueError: If the file is not a valid BAM file.
        OSError: If there is an error reading the file.
    """

    try:
        b = pysam.AlignmentFile(path, aln_format)
        header = b.header
        if 'HD' in header and 'SO' in header['HD']:
            sort_order = header['HD']['SO']
            if sort_order == 'coordinate':
                pass
                # print(f"Alignment file sort order: {sort_order}")
            else:
                print(f"Alignment file sort order: {sort_order}. It should be sorted by \'coordinate\'!!")
                print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
                sys.exit()
        else:
            print("No sort order specified in the header.")
            print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
            sys.exit()
        b.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid alignment file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()


def bam_check_tags(bam_file, args):
    """
    Check if the BAM files have the required tags for inSTRbility.

    Args:
        bam_file (str): Path to the BAM file.
        args (argparse.Namespace): Parsed command line arguments.

    Returns:
        None
    Raises:
        ValueError: If the BAM file does not contain the required tags.
    """

    reads_sampled = 0 # sample set of reads to look for tags
    aln_file = pysam.AlignmentFile(bam_file, args.aln_format)
    read_length = 0
    cs_tag = False; md_tag = False; cigar_tag = False
    for read in aln_file.fetch():
        # 0x400 - read is PCR or optical duplicate
        # 0x100 - not primary alignment
        if (read.flag & 0X400) or (read.flag & 0X100): continue
        reads_sampled += 1
        cigar = read.cigarstring
        # ?? I think the read length and also if the data is SRS or LRS should be provided by the user
        # ?? Also should read_length be summed?
        read_length += read.query_length

        if not cs_tag and read.has_tag('cs'):
            print("CS tag detected. Processing using CS tag...\n")
            cs_tag = True

        elif not cigar_tag and ((cigar!=None) and (('X' in cigar) or ('=' in cigar))):
            print("CIGAR(X/=) tag detected. Processing using CIGAR(X/=) tag...\n")
            cigar_tag = True

        elif not md_tag and read.has_tag('MD'):
            print("MD tag detected. Processing using MD tag...")
            print("Include CS tag or CIGAR tag with 'X/=' for faster processing.\n")
            md_tag = True

        if reads_sampled > 100:
            if not (cs_tag or md_tag or cigar_tag):
                print(f"No tags detected in {bam_file.split('/')[-1]}. Processing without tags...")
                print("Include the CS tag, MD tag, or CIGAR tag with 'X/=' for faster processing.\n")
            break
            # sys.exit()
    aln_file.close()

    # if average read length is less than 350bp consider it as short read data
    # ?? also this needs to be fixed as error model we might consider for short read and long read data is different
    if read_length/reads_sampled < 350:
        print('Short reads detected... Processing in short-read mode.')
        args.srs = True
    else: print('Long reads detected... Processing in long-read mode.')



def tabix_check(path):
    """
    Check if the provided regions file is valid and indexed.

    Args:
        path (str): Path to the regions file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid tabix file.
        OSError: If there is an error reading the file.
    """

    try:
        t = pysam.TabixFile(path)
        t.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid tabix file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()


def split_bedfile(tbx, contigs, total_loci, threads):
    """
    Split the regions file into chunks for parallel processing.

    Args:
        tbx (pysam.TabixFile): Tabix file object for the regions file.
        contigs (list): List of contigs to be processed.
        total_loci (int): Total number of loci in the regions file.
        threads (int): Number of threads to split the work into.

    Returns:
        list: A list of tuples, each containing a chunk of contigs and their start and
              end coordinates.
    """

    split_point = total_loci // threads
    # split_point is 0 when the total_loci is less than the number of threads
    # this is a rare case with a bed file with very few loci; all these loci will be handled by a single thread
    if split_point == 0: split_point = total_loci

    coordinate_splits = []
    line_count = 0
    current_split = []
    for contig in contigs:
        init = False
        for row in tbx.fetch(contig):
            line_count += 1
            if init == 0:
                fields = row.split('\t')
                chrom  = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                init = True

            if len(coordinate_splits) < threads-1:
                if line_count % split_point == 0:
                    end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
                    current_split.append([chrom, start_coord, end_coord])
                    coordinate_splits.append(tuple(current_split))
                    line_count = 0
                    current_split = []
                    init = False
        if init:
            end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
            current_split.append([chrom, start_coord, end_coord])
    coordinate_splits.append(tuple(current_split))

    tbx.close()

    return coordinate_splits


def main():
    """Entry point for the inSTRbility program."""

    start_time = timeit.default_timer()
    args = parse_args()

    for arg in vars(args):
        print (arg, getattr(args, arg))
    print('\n')

    # checking the input alignment format. default is bam
    if    args.aln_format == 'cram': args.aln_format = 'rc'
    elif  args.aln_format == 'sam':  args.aln_format = 'r'
    else: args.aln_format = 'rb'

    # checking the input files
    fasta_check(args.ref)
    for bam in args.bam:
        bam_check(bam, args.aln_format)

    # tabix compresses the bed file and indexes it
    pysam.tabix_compress(args.bed, f'{args.bed}.gz', force=True)
    pysam.tabix_index(f'{args.bed}.gz', seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)
    args.bed = f'{args.bed}.gz'

    tbx  = pysam.Tabixfile(args.bed)

    # the contigs to be genotyped
    total_loci = 0
    if not args.contigs:
        args.contigs = sorted(tbx.contigs)
        for row in tbx.fetch(): total_loci += 1
    else:
        args.contigs = sorted(args.contigs)
        for contig in args.contigs:
            for row in tbx.fetch(contig): total_loci += 1

    # karyotype list holds the karyotype of each sample
    if not args.karyotype:
        karyotype_list = [False]*len(args.bam)
    else:
        karyotype_list = [i=='XY' for i in args.karyotype]

    threads = args.threads
    coordinate_splits = split_bedfile(tbx, args.contigs, total_loci, args.threads)

    # each bam file is processed in separately
    for sidx, bam in enumerate(args.bam):

        print(f"Processing sample {bam.split('/')[-1]}\n")

        # checks for the necessary tags in the bam file
        bam_check_tags(bam, args)

        if threads > 1:
            thread_pool = list()
            # initializing threads
            for tidx in range(threads):
                coordinate_range = coordinate_splits[tidx]

                thread = Process(target = extract_reads, args = (args, bam, coordinate_range, tidx, karyotype_list[sidx]))
                thread.start()
                thread_pool.append(thread)

            # joining Threads
            for tidx, thread in enumerate(thread_pool):
                thread.join()
            # emptying thread_pool
            thread_pool.clear()

            out = open(f'{args.output}', 'a')
            print('Concatenating thread outputs!', file=sys.stderr)
            for tidx in range(threads)[1:]:
                thread_out = f'{args.output}_thread_{tidx}.out'
                with open(thread_out, 'r') as fh:
                    # if tidx!=0: next(fh)
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        print(*repeat_info, file=out, sep='\t')
                os.remove(thread_out)
            out.close()
            print('Concatenation completed!! ^_^', file=sys.stderr)

            if args.debug:
                log_file = open(f'{args.output}_debug.log', 'a')
                for tidx in range(threads)[1:]:
                    thread_log_out = f'{args.output}_debug_{tidx}.log'
                    with open(thread_log_out, 'r') as fh:
                        for line in fh:
                            log_info = line.strip()
                            print(log_info, file=log_file)
                    os.remove(thread_log_out)
                log_file.close()

        # inSTRbility is run in single threaded mode
        else:
            coordinate_range = coordinate_splits[0]
            # function extracts the reads from the bam file and processes them to analyse each locus
            extract_reads(args, bam, coordinate_range, 0, karyotype_list[sidx])

    time_now = timeit.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))

if __name__ == '__main__':
    main()
