#!/usr/bin/env python3
"""Step 4 - assemble Table 2 and check it against the committed baseline.

This is the acceptance step: if the numbers below match
results/table2_candidate_cpgs.tsv, the pipeline reproduced the published table.

  work/smr_input/smr_multi_finngen_pe.msmr    step 3, all SMR / HEIDI quantities
  work/smr_hits_finngen_pe.tsv                step 2, the PACE-reported mQTL effect
  -> work/table2_candidate_cpgs.tsv
  -> PASS/FAIL against results/table2_candidate_cpgs.tsv

Two details about where each number comes from:

  * SMR quantities (b_SMR, SE, P_SMR, P_SMR-multi, P_HEIDI, n SNP) are read from the
    .msmr file, not from the step-2 table. The binary writes its inputs with six
    significant digits, so its P_SMR differs from the Python screen in the third digit
    (1.57177e-07 vs 1.56420e-07 for cg06177599). Table 2 reports the binary's values,
    because those are the ones HEIDI and SMR-multi were computed with.
  * The mQTL effect and its P value are the PACE-reported values, carried through
    step 2. They are not the binary's p_eQTL, which is recomputed from the rounded
    b and se and can be off by more than 10% for very small P.

locus, tier and ld_r2_index are author-assigned quantities (LD r-squared to the
index CpG was computed separately from the 1000G reference panel); they are
constants here because they are reported values, not outputs of this pipeline.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd
import numpy as np

# reported annotation: nearest gene, cytoband, LD r2 to the index CpG, evidence tier
ANNOTATION = {
    "cg07980026": ("MTHFR", "1p36.22", 1.000, 1),
    "cg06177599": ("CLCN6", "1p36.22", 0.945, 3),
    "cg20620484": ("CLCN6", "1p36.22", 0.937, 3),
    "cg07027493": ("FLT1", "13q12.3", 1.000, 2),
    "cg06204421": ("ZNF831", "20q13.33", 1.000, 3),
}
ORDER = ["cg07980026", "cg06177599", "cg20620484", "cg07027493", "cg06204421"]
COLS = ["cpg", "nearest_gene", "locus", "top_snp", "mqtl_beta", "mqtl_p",
        "b_smr", "ci_low", "ci_high", "p_smr", "p_smr_multi", "p_heidi",
        "n_snp_heidi", "heidi_class", "ld_r2_index", "tier"]
# compared in log10 space, because the baseline stores them to 3 significant digits
PVAL = ["mqtl_p", "p_smr", "p_smr_multi"]
NUMERIC = ["mqtl_beta", "b_smr", "ci_low", "ci_high", "p_heidi", "n_snp_heidi",
           "ld_r2_index"]


def heidi_class(p: float) -> str:
    if p > 0.10:
        return "not_rejected"
    return "nominal" if p > 0.05 else "rejected"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--msmr", default="work/smr_input/smr_multi_finngen_pe.msmr")
    ap.add_argument("--hits", default="work/smr_hits_finngen_pe.tsv")
    ap.add_argument("--baseline", default="results/table2_candidate_cpgs.tsv")
    ap.add_argument("--out", default="work/table2_candidate_cpgs.tsv")
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--atol", type=float, default=5e-4)
    ap.add_argument("--plog", type=float, default=0.01,
                    help="allowed |difference| in log10 for P values")
    a = ap.parse_args()

    for path in (a.msmr, a.hits):
        if not os.path.exists(path):
            print(f"missing {path}\nrun steps 2 and 3 first", file=sys.stderr)
            return 2
    m = pd.read_csv(a.msmr, sep="\t").set_index("probeID")
    hits = pd.read_csv(a.hits, sep="\t").set_index("cpg")
    missing = set(ORDER) - set(m.index)
    if missing:
        print(f"{a.msmr} is missing CpGs: {sorted(missing)}", file=sys.stderr)
        return 1

    rows = []
    for c in ORDER:
        r, h = m.loc[c], hits.loc[c]
        gene, locus, r2, tier = ANNOTATION[c]
        rows.append({
            "cpg": c, "nearest_gene": gene, "locus": locus,
            "top_snp": r.topSNP, "mqtl_beta": h.beta, "mqtl_p": h.pace_p,
            "b_smr": r.b_SMR,
            "ci_low": r.b_SMR - 1.96 * r.se_SMR,
            "ci_high": r.b_SMR + 1.96 * r.se_SMR,
            "p_smr": r.p_SMR, "p_smr_multi": r.p_SMR_multi,
            "p_heidi": r.p_HEIDI, "n_snp_heidi": int(r.nsnp_HEIDI),
            "heidi_class": heidi_class(r.p_HEIDI),
            "ld_r2_index": r2, "tier": tier})
    t = pd.DataFrame(rows)[COLS]

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    t.to_csv(a.out, sep="\t", index=False, float_format="%.3g")
    print(t.to_string(index=False))
    print(f"\nwrote {a.out}")

    # ------------------------------------------------------------- baseline
    if not os.path.exists(a.baseline):
        print(f"\nno baseline at {a.baseline}; skipping the comparison")
        return 0
    base = pd.read_csv(a.baseline, sep="\t", encoding="utf-8-sig").set_index("cpg")
    got = t.set_index("cpg")
    if set(base.index) != set(got.index):
        print(f"\nFAIL: CpG sets differ: {sorted(set(base.index) ^ set(got.index))}")
        return 1

    worst, bad = 0.0, []
    for c in base.index:
        if got.loc[c, "heidi_class"] != base.loc[c, "heidi_class"]:
            bad.append((c, "heidi_class", base.loc[c, "heidi_class"],
                        got.loc[c, "heidi_class"], float("inf")))
        for col in PVAL:
            b, g = float(base.loc[c, col]), float(got.loc[c, col])
            d = abs(np.log10(g) - np.log10(b))
            worst = max(worst, d)
            if d > a.plog:
                bad.append((c, col, b, g, d))
        for col in NUMERIC:
            b, g = float(base.loc[c, col]), float(got.loc[c, col])
            err = abs(g - b)
            rel = err / max(abs(b), 1e-300)
            worst = max(worst, rel)
            if rel > a.rtol and err > a.atol:
                bad.append((c, col, b, g, rel))

    print(f"\n== comparison with {a.baseline} ==")
    print(f"    rows {len(base)}; cells {len(base) * (len(NUMERIC) + len(PVAL) + 1)}")
    print(f"    P values: |delta log10| <= {a.plog:g}")
    print(f"    other cells: {a.rtol:g} relative or {a.atol:g} absolute "
          f"(the baseline rounds P_HEIDI to three decimals)")
    if bad:
        print(f"    FAIL: {len(bad)} cell(s) outside tolerance")
        for c, col, b, g, e in bad:
            print(f"      {c} {col}: baseline {b} vs obtained {g} (metric {e:.2e})")
        return 1
    print("    PASS: reproduced")
    return 0


if __name__ == "__main__":
    sys.exit(main())
