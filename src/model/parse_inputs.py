from tqdm import tqdm

import sys
import gzip
import cyvcf2


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_input_tsv( input_path: str, genotypes:  dict | None) -> list:
    """
    Stream-parse the per-read TSV into a list of locus dicts.

    @param input_path the path to the TSV file (can be gzipped)
    @param genotypes  dict of {locus_key: [allele_hap0, allele_hap1]} from VCF

    Each dict: {chrom, start, end, motif, haplotypes: {hap: [lengths_ru]},
                genotypes: array|None}
    Lengths are converted to repeat units (bp / motif_size).
    """
    ins = (gzip.open(input_path, "rt") if input_path.endswith(".gz") else open(input_path, "rt"))

    data:  list = []
    info:  dict = {}
    haps:  dict = {}
    prev_key: str | None = None

    pbar = tqdm(unit="rows", unit_scale=True, ncols=80, smoothing=0.1, position=1, desc="Reading TSV")

    def _flush(data, info, key, haps, genotypes):
        """
        Flush the current locus data into the output list.

        @param data      the output list of locus dicts
        @param info      dict of {chrom, start, end, motif} for the current locus
        @param key       the current locus key
        @param haps      dict of {hap: [lengths_ru]} for the current locus
        @param genotypes dict of {locus_key: [allele_hap0, allele_hap1]} from VCF
        """
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

        fields = line.rstrip("\n").split("\t")
        (chrom, start, end, motif, read_id, haplotype, length_bp, allele, avg_meth, meth_bases) = fields

        start         = int(start)
        end           = int(end)
        ml            = len(motif)
        ref_length_bp = (end - start)
        pbar.update(1)

        key = f"{chrom}:{start}-{end}_{motif}"
        info.setdefault(key, {"chrom": chrom, "start": start, "end": end, "motif": motif})

        hap   = int(haplotype)
        units = float(length_bp) / ml

        if key != prev_key and prev_key is not None:
            _flush(data, info, prev_key, haps, genotypes)
            haps = {}
            del info[prev_key]

        haps.setdefault(hap, []).append(units)
        prev_key = key

    if haps and prev_key is not None:
        _flush(data, info, prev_key, haps, genotypes)

    pbar.close()
    ins.close()

    print(f"Loaded data for {len(data)} loci from {input_path}", file=sys.stderr)
    return data


def load_genotypes_from_vcf(vcf_path: str) -> dict:
    """
    Load genotypes from a VCF file into a dict of {locus_key: [founder_length_hap0, founder_length_hap1]}.
    Locus key is of the form "chrom:start-end_motif".
    
    @param vcf_path the path to the VCF file
    @return dict of {locus_key: [founder_length_hap0, founder_length_hap1]}
    """

    pbar = tqdm(unit="rows", unit_scale=True, ncols=80, smoothing=0.1, position=1, desc="Reading VCF")

    genotypes = {}
    for v in cyvcf2.VCF(vcf_path):
        if v.FILTER not in ("PASS", None):
            continue

        motif = v.INFO.get("MOTIF", "")
        if not motif: continue

        al = v.format("AL")
        if al is not None and len(al):
            key = f"{v.CHROM}:{v.start}-{v.INFO.get('END')}_{motif}"
            genotypes[key] = al[0] / len(motif)
        pbar.update(1)

    pbar.close()
    print(f"Loaded {len(genotypes)} VCF genotypes from {vcf_path}", file=sys.stderr)
    return genotypes