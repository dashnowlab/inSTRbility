#!/usr/bin/env python3

from nbgeom_modelling import (
    analyse_genome_wide, write_tsv, Regime1GlobalParams
)

import time
import argparse
import sys
import gzip
from tqdm import tqdm
from pathlib import Path

import cyvcf2
import numpy as np


def _parse_args():
    p = argparse.ArgumentParser(
        description="NB-Geometric somatic instability analysis"
    )
    p.add_argument("-i", "--input", required=True,   help="Input tsv file with read allele lengths for catalog loci")
    p.add_argument("-v", "--vcf",   required=False,  help="VCF file with genotype information for each locus."\
                                                          "If not provided, modal length for each locus will be used.")
    p.add_argument("-o", "--output", help="Output tsv file for global params", default=None)
    return p.parse_args()


def load_genotypes_from_vcf(vcf_path):
    """
    Load genotypes from a VCF file into a dictionary.
    The keys are formatted as "CHROM:start-END_MOTIF" and the values are the corresponding genotypes.
    @param vcf_path: Path to the VCF file.
    @return: Dictionary of genotypes.
    """
    start_time = time.time()
    genotypes = {}
    cyvcf_reader = cyvcf2.VCF(vcf_path)
    for variant in cyvcf_reader:
        if variant.FILTER not in ("PASS", None):
            continue
        key = f"{variant.CHROM}:{variant.start}-{variant.INFO.get('END')}_{variant.INFO.get('MOTIF')}"
        genotypes[key] = variant.format("AL")[0]

    print(f"Loaded genotypes for {len(genotypes)} loci in {time.time() - start_time:.2f} seconds", file=sys.stderr)
    return genotypes


if __name__ == "__main__":
    args = _parse_args()
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None

    # Global params: either from Pass 1 estimation or defaults
    gp = Regime1GlobalParams.from_stratum("str_short")   # quick start
    # gp = estimate_regime1_global(reference_loci)        # proper Pass 1

    if args.output is None:
        args.output = Path(args.input).parent / (Path(args.input).stem + "_instability.tsv")
    fout = open(args.output, 'wt')

    header = ["bam_id", "chrom", "start", "end", "motif", "haplotype", "founder_length", "length_in_units",
              "avg", "n_reads", "r_e", "r_c", "p_plus", "p_minus", "mu_expansion",
              "mu_contraction", "instability_index", "net_bias", "p_expansion",
              "expansion_rate", "contraction_rate", "net_mutation_rate", "p_any_mutation"]
    print("\t".join(header), file=fout)

    progress = tqdm(unit="loci", unit_scale=True, ncols=80, smoothing=0.1, desc="Processing loci")
    if genotypes is not None:
        progress.total = len(genotypes)

    haplotype_alens = {}
    haplotypes = []
    prev_key = None
    fh = gzip.open(args.input, 'rt') if args.input.endswith('.gz') else open(args.input, 'rt')

    for line in fh:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        chrom, start, end, motif, read_id, haplotype, length, allele, avg_meth, meth_bases = line
        locus_key = f'{chrom}:{start}-{end}_{motif}'
        haplotype = int(haplotype)
        length = float(length)
        units  = length / len(motif)
        haplotypes = {}
        if haplotype not in haplotypes:
            haplotypes[haplotype] = []
        haplotypes[haplotype].append(float(units))

        if locus_key != prev_key and prev_key is not None:
            if haplotypes:
                if samples[prev_key].startswith("FXPM5006"):
                    alens = []
                    for h, hap in enumerate(sorted(haplotypes)):
                        alens += haplotype_alens[hap]
                    fl = _modal_length(np.array(alens))
                    ru = fl
                    avg = average_median_delta(alens)
                    print(samples.get(prev_key, "unknown"), min(alens), max(alens), avg, fl, ru, alens)
                    process_locus(chrom, start, end, motif, h, alens, fl, ru, out, samples.get(prev_key, "unknown"), avg)
                else:
                    for h, hap in enumerate(sorted(haplotypes)):
                        alens = haplotype_alens[hap]
                        fl = genotypes[prev_key][h] if prev_key in genotypes else np.median(alens)
                        ru = fl
                        avg = average_median_delta(alens)
                        # ru = float(fl)/len(motif)
                        # if ru >= 3 and fl >= 20:
                        process_locus(chrom, start, end, motif, h, alens, fl, ru, out, samples.get(prev_key, "unknown"), avg)
                progress.update(1)
            haplotype_alens = {}
            haplotypes      = []
        prev_key = key
    if haplotypes:
        for h, hap in enumerate(sorted(haplotypes)):
            alens = haplotype_alens[hap]
            fl = genotypes[prev_key][h] if prev_key in genotypes else np.median(alens)
            avg = average_median_delta(alens)
            process_locus(chrom, start, end, motif, h, alens, fl, fl, out, samples.get(prev_key, "unknown"), avg)
    fh.close()
    if out != sys.stdout: out.close()
    progress.close()

    results = analyse_genome_wide(
        catalog,
        global_params=gp,
        founder_lengths=founder_lengths,   # optional but recommended
        mitotic_generations=60.0,          # adult blood; omit if unknown
        n_ppp=300,                         # PPP simulations per locus
        verbose=True,
    )
    # Load the catalog from the input TSV file
    catalog = load_catalog_from_tsv(args.input, genotypes)
    # Perform the analysis and write results
    results = analyse_genome_wide(catalog, global_params=gp, founder_lengths=founder_lengths)
    write_tsv(results, args.output or "instability_results.tsv", include_per_div=True)
    write_tsv(results, "instability_results.tsv", include_per_div=True)