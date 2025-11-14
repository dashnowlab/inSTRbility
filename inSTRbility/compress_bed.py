import pysam
from filelock import FileLock
import sys
import os
import errno


def safe_tabix_compress(in_bed, out_bgz, force=False):
    # Open output file for writing
    lock = FileLock(in_bed + ".lock")
    with lock:
        # tabix compresses the bed file and indexes it
        try: pysam.tabix_compress(in_bed, out_bgz, force=False)
        except OSError: pass
        print(f'Compressed {in_bed} to {out_bgz}')
        try: pysam.tabix_index(out_bgz, seq_col=0, start_col=1, end_col=2, zerobased=True, force=False)
        except OSError: pass
    os.remove(in_bed + ".lock")


if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[1] + ".gz"
    safe_tabix_compress(input, output, force=True)
