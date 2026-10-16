#!/usr/bin/env python3
"""Validate the UMI-collision (saturation) correction on a real library, using ESAT's per-key output.

Input: <exp>.umi.distributions.txt[.gz] from an ESAT run with -scPrep -umiMode position (one row per key:
barcode, gene, UMI@5'position, reads). Two tests, both with real position profiles, UMI frequencies and PCR
duplication; none (keys, no correction) and exact (ESAT -umiSaturation exact):

  truncation   genes nearly unsaturated with the full UMIs (k/K < 0.2) give the truth; every UMI is cut to
               5, 4, 3 nt and exact must recover that truth (the validation of Agyemang, Irizarry & Baharav,
               J Comput Biol 2026, on shortened UMIs).
  robustness   exact (and keys) against its own full-UMI estimate, for all genes: how much does the answer change when
               UMIs saturate 4x, 16x, 64x more easily? R^2 around the identity (log2 scale), bias and spread.

Usage: validate_saturation.py <umi.distributions.txt[.gz]> <out_prefix> [lengths, default 5,4,3]
Writes <out_prefix>.truncation.tsv, <out_prefix>.robustness.tsv (per gene) and prints the summaries.
"""
import sys
import numpy as np
import pandas as pd
import umi_saturation as us

path, out = sys.argv[1], sys.argv[2]
lengths = tuple(int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "5,4,3").split(","))
keys = us.read_keys(path)
pd.set_option("display.width", 250)

tr = us.truncation_test(keys, lengths=lengths, truth_max_frac=0.2)
tr.to_csv(out + ".truncation.tsv", sep="\t")
tr["k/K"] = pd.cut(tr["k_over_K"], [0, 0.2, 0.5, 0.8, 1.5, np.inf]).astype(str)
rows = []
for (L, b), d in tr.groupby(["L", "k/K"]):
    r = {"UMI length": L, "k/K": b, "genes": len(d)}
    for m, lab in [("k", "none"), ("exact", "exact")]:
        e = np.log2(d[m] / d["truth"])
        r[f"{lab} bias"], r[f"{lab} sd"] = round(float(e.median()), 3), round(float(1.4826 * (e - e.median()).abs().median()), 3)
    rows.append(r)
print("== truncation test: log2(estimate / truth), median (bias) and robust SD")
print(pd.DataFrame(rows).to_string(index=False))

sc = us.self_consistency(keys, lengths=lengths)
sc.to_csv(out + ".robustness.tsv", sep="\t")
print("\n== robustness: estimate with shortened UMIs vs the same method with the full UMIs")
print(us.self_consistency_summary(sc).to_string(index=False))
