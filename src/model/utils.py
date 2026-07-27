# !/usr/bin/env python3

stratum_labels = {1: "1bp", 2: "2bp", 3: "3bp", 4: "4bp", 5: "5bp", 6: "6bp", 7: "7+bp", 20: ">20bp"}
def get_stratum(motif_len: int) -> int:
    """
    Motif length stratum for global parameter estimation. Returns 3, 6, or 7.

    @param motif_len motif length in bp
    @return stratum int
    """
    if motif_len <= 6:   return motif_len
    if motif_len <= 20:  return 7
    else:                return 20
