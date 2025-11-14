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