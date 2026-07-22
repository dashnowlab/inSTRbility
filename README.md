# inSTRbility


inSTRbility is a toolkit to analyse somatic instability at tandem repeat loci from whole genome sequencing datasets.

The basic pipeline of the tool includes fetching reads mapping to repeat region and recording variations within the 
repeat region along with SNPs outside of the repeat. The variations falling within the repeat region contribute towards
calculation of the allele length in the read. Based on the SNP the reads are segregated into two haplogroups. The allele
length for each group is calculated as the mode of the group. For each read instability is calculated as the mean absolute
deviation from the allele length.

<b>NOTE:</b> The tool currently works on long read sequencing datasets including PacBio and ONT. 


## Usage

### Printing help
```bash
$ python ./inSTRbility/core.py -h
```

### Basic usage
```bash
$ python ./inSTRbility/core.py -ref [fasta] -bed [regions_file] -bam [aln_file] -o [output_file] --reads-out
```

### Scripts for extracting reads
- __init__.py
- core.py
- cigar_utils.py
- cstag_utils.py
- md_utils.py
- genotype_utils.py
- cigar_utils.py
- phasing_utils.py
- instabillity_utils.py
- process_reads.py
- locus_utils.py
- operation_utils.py
- version.py

### Scripts for modelling instability
- calc_instability.py
- parse_inputs.py
- suggest_re_rc_grid.py
- nbgeom_modelling.py
- compress_bed.py
- calibrate.py
- check_lock.py
- multisample-analysis.py
