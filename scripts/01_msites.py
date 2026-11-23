#!/usr/bin/env python3
"""Step 1 - collapse the public PACE top-SNP export to one instrument per CpG.

The public PACE placental cis-mQTL release is distributed as one row per
significant SNP-CpG pair. Every downstream analysis in the paper works on the
CpG level, so this step keeps, for each CpG, the single pair with the strongest
association and reports how many rows each filter removes.

Input   data/raw/pace_cis_mqtl_top_snps.tsv   (19-column PACE export, manual)
Output  work/msites.tsv                       288,760 rows

The 288,760 figure is the "CpG-level mSites" quoted in Methods and Fig. 2c. It is
not the number of released pairs: the export contains 655,543 of those.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
from scipy import special

SNP_RE = r"(\d{1,2}|X|Y):(\d+):([ACGT]+)_([ACGT]+)"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/raw/pace_cis_mqtl_top_snps.tsv")
    ap.add_argument("--out", default="work/msites.tsv")
    ap.add_argument("--expect-msites", type=int, default=288760)
    ap.add_argument("--no-check", action="store_true",
                    help="report the counts but do not enforce them")
    a = ap.parse_args()

    pace = pd.read_csv(a.input, sep="\t", dtype=str, low_memory=False)
    n_raw = len(pace)

    # ---- filter 1: the snp column must be chr:pos:A1_A2 -------------------
    ok = pace["snp"].astype(str).str.fullmatch(SNP_RE)
    n_bad = int((~ok).sum())
    pace = pace[ok].copy()

    # ---- filter 2: exact duplicate cpg-snp rows ---------------------------
    n_dup = int(pace.duplicated(["cpg", "snp"]).sum())
    pace = pace[~pace.duplicated(["cpg", "snp"])].copy()

    sp = pace["snp"].str.split(":", expand=True)
    pace["snp_chr"] = sp[0]
    pace["snp_pos"] = sp[1].astype(np.int64)
    # the export stores the effect allele under reference_allele and its
    # frequency under eaf; both are verified against GWAS allele frequencies
    # downstream in step 2
    for src, dst in [("beta", "beta"), ("se", "se"), ("z", "z"),
                     ("p-value", "p"), ("n_studies", "n_studies")]:
        pace[dst] = pd.to_numeric(pace[src], errors="coerce")
    n_nan = int(pace[["beta", "z", "n_studies"]].isna().any(axis=1).sum())
    pace = pace.dropna(subset=["beta", "z", "n_studies"]).copy()
    pace["n_studies"] = pace["n_studies"].astype(int)

    # -log10(P) recomputed from z so that it stays finite when the reported
    # p-value underflows to zero
    pace["neg_log10p"] = -(np.log(2.0) + special.log_ndtr(-pace["z"].abs())) / np.log(10.0)

    # ---- collapse to the strongest pair per CpG ---------------------------
    msites = (pace.sort_values("neg_log10p", ascending=False)
                  .drop_duplicates("cpg")
                  .reset_index(drop=True))

    cols = ["cpg", "snp", "snp_chr", "snp_pos", "reference_allele",
            "other_allele", "eaf", "beta", "se", "z", "p", "neg_log10p",
            "n_studies"]
    msites = msites[cols]
    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    msites.to_csv(a.out, sep="\t", index=False, float_format="%.10g")

    print(f"input rows                  {n_raw:>10,}")
    print(f"  malformed snp             {-n_bad:>10,}")
    print(f"  duplicate cpg-snp         {-n_dup:>10,}")
    print(f"  unparseable numeric       {-n_nan:>10,}")
    print(f"retained pairs              {n_raw - n_bad - n_dup - n_nan:>10,}")
    print(f"CpG-level mSites            {len(msites):>10,}   -> {a.out}")

    if not a.no_check and len(msites) != a.expect_msites:
        print(f"\nFAIL: expected {a.expect_msites:,} CpG-level mSites, got {len(msites):,}.\n"
              "The PACE export is a live resource; if you downloaded it on a different\n"
              "date the count can differ by a handful. Re-check data/README.md, or pass\n"
              "--no-check if you intentionally want to continue.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
