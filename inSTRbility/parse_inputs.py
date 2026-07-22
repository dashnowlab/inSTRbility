from tqdm import tqdm

import sys
import time
import gzip
import cyvcf2


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_input_tsv( input_path: str, genotypes:  dict | None) -> list:
    """
    Stream-parse the per-read TSV into a list of locus dicts.

    @param input_path the path to the TSV file (can be gzipped)
    @param genotypes  dict of {locus_key: [founder_length_hap0, founder_length_hap1]} from VCF

    Each dict: {chrom, start, end, motif, haplotypes: {hap: [lengths_ru]},
                genotypes: array|None}
    Lengths are converted to repeat units (bp / motif_size).
    """
    ins = (gzip.open(input_path, "rt")
           if input_path.endswith(".gz")
           else open(input_path, "rt"))

    data:  list = []
    info:  dict = {}
    haps:  dict = {}
    prev_key: str | None = None

    pbar = tqdm(unit="rows", unit_scale=True, ncols=80, smoothing=0.1, position=1, desc="Reading TSV")

    def _flush(data, info, key, haps, genotypes):
        data.append({
            "chrom":      info[key]["chrom"],
            "start":      info[key]["start"],
            "end":        info[key]["end"],
            "motif":      info[key]["motif"],
            "haplotypes": dict(haps),
            "genotypes":  (genotypes.get(key) if genotypes else None),
        })

    for line in ins:
        if line.startswith("#"): continue

        parts = line.rstrip("\n").split("\t")
        (chrom, start, end, motif, _read_id, haplotype, length_bp, _allele, _avg_meth, _meth_bases) = parts

        start  = int(start)
        end    = int(end)
        ml     = len(motif)
        span   = end - start
        pbar.update(1)

        key = f"{chrom}:{start}-{end}_{motif}"
        info.setdefault(key, {"chrom": chrom, "start": start, "end": end, "motif": motif})

        hap  = int(haplotype)
        unit = float(length_bp) / ml

        if key != prev_key and prev_key is not None:
            _flush(data, info, prev_key, haps, genotypes)
            haps = {}
            del info[prev_key]

        haps.setdefault(hap, []).append(unit)
        prev_key = key

    if haps and prev_key is not None:
        _flush(data, info, prev_key, haps, genotypes)

    pbar.close()
    ins.close()
    return data

def load_genotypes_from_vcf(vcf_path: str) -> dict:
    """
    Load genotypes from a VCF file into a dict of {locus_key: [founder_length_hap0, founder_length_hap1]}.
    Locus key is of the form "chrom:start-end_motif".
    
    @param vcf_path the path to the VCF file
    @return dict of {locus_key: [founder_length_hap0, founder_length_hap1]}
    """
    t0, genotypes = time.time(), {}
    for v in cyvcf2.VCF(vcf_path):
        if v.FILTER not in ("PASS", None):
            continue
        motif = v.INFO.get("MOTIF", "")
        if not motif: continue
        al = v.format("AL")
        if al is not None and len(al):
            key = f"{v.CHROM}:{v.start}-{v.INFO.get('END')}_{motif}"
            genotypes[key] = al[0] / len(motif)
    print(f"Loaded {len(genotypes)} VCF genotypes in {time.time()-t0:.1f}s",
          file=sys.stderr)
    return genotypes