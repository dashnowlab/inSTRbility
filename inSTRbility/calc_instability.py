#!/usr/bin/env python3

from xml.parsers.expat import errors

from nbgeom_modelling import analyse_haplotype, _modal_length, Regime1GlobalParams, estimate_regime1_global

import time
import argparse
import sys
import gzip
from tqdm import tqdm
from pathlib import Path
from multiprocessing import Pool

import cyvcf2
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed


def _parse_args():
    p = argparse.ArgumentParser(
        description="NB-Geometric somatic instability analysis"
    )
    p.add_argument("-i", "--input", required=True,   help="Input tsv file with read allele lengths for catalog loci")
    p.add_argument("-v", "--vcf",   required=False,  help="VCF file with genotype information for each locus."\
                                                          "If not provided, modal length for each locus will be used.")
    p.add_argument("-o", "--output", help="Output tsv file for global params", default=None)
    p.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for parallel processing")
    return p.parse_args()


def average_median_delta(alens, founder_length=None):
    """
    Calculate the average median delta for a list of allele lengths.
    The average median delta is defined as the mean of the absolute differences between each allele length and the median allele length.
    @param alens: List of allele lengths.
    @param founder_length: The founder length for the locus.
    @return: Average median delta.
    """
    lengths = np.asarray(alens)
    if founder_length is not None: founder_length = _modal_length(alens)
    deltas = (np.abs(lengths - founder_length))
    return np.mean(deltas)


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
        genotypes[key] = variant.format("AL")[0]/len(variant.INFO.get('MOTIF'))

    print(f"Loaded genotypes for {len(genotypes)} loci in {time.time() - start_time:.2f} seconds", file=sys.stderr)
    return genotypes


def _worker_function(chunk_data, fout, t):
    """
    Worker function to process a single locus and write results to the output file.
    @param locus_data: Dictionary containing locus information and haplotype data.
    @param fout: Output file handle.
    """
    header = ["chrom", "start", "end", "motif", "haplotype", "founder_length", "length_in_units",
              "avg_deviation", "n_reads", "r_e", "r_c", "p_plus", "p_minus", "mu_expansion", "mu_contraction","instability_index",
              "net_bias", "p_expansion", "expansion_rate", "contraction_rate", "net_mutation_rate", "p_any_mutation", 
              "ks_statistic", "ad_statistic", "ppp_variance", "ppp_skewness", "observed_var", "fitted_var", "dispersion_ratio",
              "chisq_statistic", "chisq_pvalue", "chisq_df", "chisq_n_bins", "regime2_recommended","fit_quality", "fit_quality_adjusted"]

    # t_progress = tqdm(unit="loci", total=len(chunk_data), unit_scale=True, ncols=80,
    #                   smoothing=0.1, desc="Processing loci (worker)", position=2*(t + 1))

    out = open(fout, 'wt')
    if t == 0:
        print("\t".join(header), file=out)

    for locus_data in chunk_data:
        chrom = locus_data["chrom"]
        start = locus_data["start"]
        end = locus_data["end"]
        motif = locus_data["motif"]
        haplotypes = locus_data["haplotypes"]
        genotypes = locus_data["genotypes"]
        for h, hap in enumerate(sorted(list(haplotypes.keys()))):
            alens = haplotypes[hap]
            result = analyse_haplotype(
                lengths = alens,
                founder_length = genotypes[h],
                repeat_unit_length=len(motif)
            )
            avg = average_median_delta(alens, founder_length=result["founder_length"])
            print(chrom, start, end, motif,
                hap, result["founder_length"], round(result["founder_length"] / len(motif), 2), avg,
                result["n_reads"], result["r_e"], result["r_c"], result["p_plus"], result["p_minus"], result["mu_N_expansion"],
                result["mu_N_contraction"], result["instability_index"], result["net_bias"], result["p_expansion"],
                result["expansion_rate"], result["contraction_rate"], result["net_mutation_rate"], result["p_any_mutation"],
                result["gof_ks_statistic"], result["gof_ad_statistic"], result["gof_ppp_variance"], result["gof_ppp_skewness"],
                result["gof_observed_var"], result["gof_fitted_var"], result["gof_dispersion_ratio"],
                result["gof_chisq_statistic"], result["gof_chisq_pvalue"], result["gof_chisq_df"], result["gof_chisq_n_bins"], result["gof_regime2_recommended"],
                result["gof_fit_quality"], result["gof_fit_quality_adjusted"],
                sep="\t", file=out)
        # t_progress.update(1)
    out.close()
    # t_progress.close()


def _multiprocess(data, fout, threads, progress):
    """
    Process the data in parallel using multiple threads.
    @param data: List of data to process.
    @param fout: Output file handle.
    @param threads: Number of threads to use.
    @param progress: tqdm progress bar object.
    """

    chunks = np.array_split(data, threads)

    fouts = [fout] + [f"{fout}.part{i}" for i in range(1, threads)]

    threads_list = []
    def on_success(_): progress.update()

    def on_error(e):
        progress.update()
        print(f'Thread failed: {e}', file=sys.stderr)
    with Pool(processes=args.threads) as pool:
        results = []
        for thread_idx in range(args.threads):
            r = pool.apply_async(
                _worker_function,
                args     = (chunks[thread_idx], fouts[thread_idx], thread_idx),
                callback = on_success,
                error_callback = on_error    # catch thread failures
            )
        pool.close()
        pool.join()
    
    if threads > 1:
        out = open(fout, 'at')
        for i in range(1, threads):
            with open(f"{fout}.part{i}", 'rt') as part_file:
                for line in part_file:
                    out.write(line)
            Path(f"{fout}.part{i}").unlink()
        out.close()

    
if __name__ == "__main__":
    args = _parse_args()
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None

    # Global params: either from Pass 1 estimation or defaults
    gp = Regime1GlobalParams.from_stratum("str_short")   # quick start
    # gp = estimate_regime1_global(reference_loci)        # proper Pass 1

    if args.output is None:
        args.output = Path(args.input).parent / (Path(args.input).stem + "_instability.tsv")

    if args.threads == 1:
        main_pbar = tqdm(unit           ="loci",
                         unit_scale     = True,
                         ncols          = 80,
                         smoothing      = 0.1,
                         position       = 3,
                         desc           = "Processing loci")
        if genotypes is not None: main_pbar.total = len(genotypes)
    else:
        main_pbar = tqdm(unit           = "chunks",
                         unit_scale     = True,
                         ncols          = 80,
                         smoothing      = 0.1,
                         total          = args.threads,
                         position       = 3,
                         desc           = "Processing chunks")

    inread_pbar = tqdm(unit           = "rows",
                       unit_scale     = True,
                       ncols          = 80,
                       smoothing      = 0.1,
                       position       = 1,
                       desc           = "Reading input")
    
    ins_fh = gzip.open(args.input, 'rt') if args.input.endswith('.gz') else open(args.input, 'rt')

    data       = []
    prev_key   = None
    haplotypes = {}
    info       = {}
    for line in ins_fh:
        if line.startswith("#"): continue
        
        line = line.strip().split("\t")
        chrom, start, end, motif, read_id, haplotype, length, allele, avg_meth, meth_bases = line
        start         = int(start)
        end           = int(end)
        motif_length  = len(motif)
        repeat_length = end - start
        inread_pbar.update(1)

        if motif_length > 6 or repeat_length < 20 or repeat_length > 500:  continue
        locus_key = f'{chrom}:{start}-{end}_{motif}'
        if locus_key not in info:
            info[locus_key] = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "motif": motif
            }
        haplotype = int(haplotype)
        length    = float(length)
        units     = length / len(motif)
        if haplotype not in haplotypes:
            haplotypes[haplotype] = []
        haplotypes[haplotype].append(float(units))

        if locus_key != prev_key and prev_key is not None:
            if haplotypes:
                data.append({
                    "chrom": info[prev_key]["chrom"],
                    "start": info[prev_key]["start"],
                    "end": info[prev_key]["end"],
                    "motif": info[prev_key]["motif"],
                    "haplotypes": haplotypes,
                    "genotypes": genotypes[prev_key] if prev_key in genotypes else None
                })
            haplotypes      = {}
            if prev_key in info: del info[prev_key]
        prev_key = locus_key
    if haplotypes:
        data.append({
            "chrom": info[prev_key]["chrom"],
            "start": info[prev_key]["start"],
            "end": info[prev_key]["end"],
            "motif": info[prev_key]["motif"],
            "haplotypes": haplotypes,
            "genotypes": genotypes[prev_key] if prev_key in genotypes else None
        })
    inread_pbar.close()

    print(f"\nLoaded {len(data)} loci for analysis", file=sys.stderr)
    
    _multiprocess(data, args.output, args.threads, main_pbar)

    ins_fh.close()
    main_pbar.close()
