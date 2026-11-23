#!/usr/bin/env python3
"""Step 2 - genome-wide single-SNP SMR screen against three outcome GWASs.

Pipeline (identical to the analysis reported in the paper):

  work/msites.tsv            288,760 CpG-level mSites (hg19)
    autosomes, single-base alleles, finite b/se/eaf
    liftOver hg19 -> hg38                      -> 288,714 mapped (46 unmapped)
    drop A/T and C/G palindromic instruments   -> 248,009 usable instruments
    inner join to GWAS on (chrom, pos38); orient alleles
    drop palindromic sites that survived the join
    one row per CpG                            -> 242,924 testable (FinnGen PE)
    b_SMR = b_GWAS / b_mQTL, delta-method SE, Bonferroni 0.05 / 248,009

Inputs   work/msites.tsv
         data/raw/finngen_R13_O15_PREECLAMPS.gz
         data/raw/finngen_R13_O15_PREEC_OR_FETGRO.gz
         data/raw/GCST90301704.h.tsv.gz
         data/ref/hg19ToHg38.over.chain
Outputs  work/msites_hg38.tsv, work/smr_<outcome>.tsv,
         work/smr_hits_<outcome>.tsv, work/stats.json

Notes
  * Harmonisation is by GRCh38 chromosome and base-pair position, not by rsID.
    The PACE export carries rsIDs as part of the composite "snp" key only, and
    rsID-based joins were used at the HEIDI / SMR-multi step (step 3) where a
    variant list had to be handed to the SMR binary.
  * The ratio estimator implemented here is the SMR estimator; it is applied in
    Python rather than by the binary so that all 242,924 pairs can be screened
    in one pass. Step 3 uses the binary for BESD construction and HEIDI.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd
from scipy import stats

PALINDROMIC = {"AT", "TA", "CG", "GC"}
AUTOSOMES = [str(i) for i in range(1, 23)]

# name, path, source format, cases, controls, label used in the paper
OUTCOMES = [
    ("finngen_pe", "finngen_R13_O15_PREECLAMPS.gz", "finngen", 9151, 258775),
    ("finngen_negctrl", "finngen_R13_O15_PREEC_OR_FETGRO.gz", "finngen", 14616, 267293),
    ("gcst90301704", "GCST90301704.h.tsv.gz", "gcst", 1728, 192399),
]

EXPECTED = {
    "msites_total": 288760,
    "lifted_hg38": 288714,
    "lift_unmapped": 46,
    "instruments_tested": 248009,
    "tested:finngen_pe": 242924,
}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def read_gwas(path: str, kind: str) -> tuple[pd.DataFrame, dict]:
    """Stream a GWAS file and return the instrument-position matches.

    Both source formats are normalised to the same canonical columns so that a
    single allele-orientation code path serves all three outcomes:

        chrom (str), pos38 (int64), gwas_eff, gwas_oth,
        gwas_beta, gwas_se, gwas_p, af_eff, rsid, nearest_gene
    """
    if kind == "finngen":
        usecols = {"#chrom": str, "pos": np.int64, "ref": str, "alt": str,
                   "rsids": str, "nearest_genes": str, "pval": np.float64,
                   "beta": np.float64, "sebeta": np.float64,
                   "af_alt": np.float64}
        reader = pd.read_csv(path, sep="\t", usecols=list(usecols), dtype=usecols,
                             chunksize=2_000_000, low_memory=False)
        renames = {"#chrom": "chrom", "pos": "pos38", "alt": "gwas_eff",
                   "ref": "gwas_oth", "beta": "gwas_beta", "sebeta": "gwas_se",
                   "pval": "gwas_p", "af_alt": "af_eff", "rsids": "rsid",
                   "nearest_genes": "nearest_gene"}
    elif kind == "gcst":
        usecols = {"chromosome": str, "base_pair_location": np.int64,
                   "effect_allele": str, "other_allele": str, "beta": np.float64,
                   "standard_error": np.float64, "p_value": np.float64,
                   "effect_allele_frequency": np.float64, "rsid": str}
        reader = pd.read_csv(path, sep="\t", usecols=list(usecols), dtype=usecols,
                             chunksize=2_000_000, low_memory=False)
        renames = {"chromosome": "chrom", "base_pair_location": "pos38",
                   "effect_allele": "gwas_eff", "other_allele": "gwas_oth",
                   "beta": "gwas_beta", "standard_error": "gwas_se",
                   "p_value": "gwas_p",
                   "effect_allele_frequency": "af_eff"}
    else:
        raise ValueError(kind)

    keep, n_rows, n_aut = [], 0, 0
    for chunk in reader:
        n_rows += len(chunk)
        chunk = chunk[chunk["chrom" if kind == "gcst" else "#chrom"].isin(AUTOSOMES)]
        if chunk.empty:
            continue
        n_aut += len(chunk)
        ch = chunk.rename(columns=renames)
        ch["chrom"] = ch["chrom"].astype(str)
        if "nearest_gene" not in ch:
            ch["nearest_gene"] = ""
        if "rsid" not in ch:
            ch["rsid"] = ""
        keep.append(ch)

    g = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    stats = {"variants_read": int(n_rows), "autosomal_variants": int(n_aut)}
    return g, stats


def smr_estimate(d: pd.DataFrame, n_tests: int) -> pd.DataFrame:
    """Ratio estimator b_SMR = b_GWAS / b_mQTL with delta-method SE."""
    d = d.copy()
    d["b_smr"] = d.gwas_beta / d.beta
    d["se_smr"] = np.sqrt(d.gwas_se ** 2 / d.beta ** 2
                          + d.gwas_beta ** 2 * d.se ** 2 / d.beta ** 4)
    d["z_smr"] = d.b_smr / d.se_smr
    d["p_smr"] = 2 * stats.norm.sf(np.abs(d.z_smr))
    d["bonferroni_thr"] = 0.05 / n_tests
    d["smr_hit"] = d.p_smr < d.bonferroni_thr
    return d.sort_values("p_smr").reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--msites", default="work/msites.tsv")
    ap.add_argument("--raw", default="data/raw")
    ap.add_argument("--chain", default="data/ref/hg19ToHg38.over.chain")
    ap.add_argument("--out", default="work")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if a headline count differs from the paper")
    a = ap.parse_args()

    try:
        from pyliftover import LiftOver
    except ImportError:
        log("pyliftover is required:  pip install -r environment/requirements.txt")
        return 2

    os.makedirs(a.out, exist_ok=True)
    stats: dict = {}

    # ---------------------------------------------------------- instruments
    log("== 1. PACE instruments: hg19 -> hg38 ==")
    ms = pd.read_csv(a.msites, sep="\t",
                     usecols=["cpg", "snp", "snp_chr", "snp_pos",
                              "reference_allele", "other_allele", "eaf",
                              "beta", "se", "z", "p", "neg_log10p", "n_studies"])
    n0 = len(ms)
    ms = ms[ms.snp_chr.astype(str).str.fullmatch(r"\d{1,2}")].copy()
    ms["chr"] = ms.snp_chr.astype(int)
    ms = ms[ms.chr.between(1, 22)].copy()
    ms = ms[ms.reference_allele.str.len().eq(1) & ms.other_allele.str.len().eq(1)]
    ms = ms.dropna(subset=["beta", "se", "eaf"])
    ms = ms[(ms.se > 0) & np.isfinite(ms.beta)]
    n_attempt = len(ms)
    log(f"    CpG-level mSites {n0:,}; autosomal instruments with finite b/se {n_attempt:,}")

    if not os.path.exists(a.chain):
        if os.path.exists(a.chain + ".gz"):
            with gzip.open(a.chain + ".gz", "rb") as f, open(a.chain, "wb") as o:
                shutil.copyfileobj(f, o)
        else:
            log(f"    missing liftOver chain: {a.chain}\n"
                "    run scripts/00_fetch_data.sh, or download hg19ToHg38.over.chain.gz\n"
                "    from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/")
            return 2

    lo = LiftOver(a.chain)
    pos38 = []
    for c, p in zip(ms.chr.to_numpy(), ms.snp_pos.to_numpy()):
        hit = lo.convert_coordinate(f"chr{c}", int(p), "+")
        pos38.append(hit[0][1] if hit else -1)
    ms["pos38"] = np.asarray(pos38, dtype=np.int64)
    n_unmapped = int((ms.pos38 < 0).sum())
    ms = ms[ms.pos38 > 0].copy()
    ms["chrom"] = ms.chr.astype(str)
    log(f"    liftOver: mapped {len(ms):,} / {n_attempt:,} "
        f"({100 * len(ms) / n_attempt:.2f}%), unmapped {n_unmapped:,}")

    pal = ms.other_allele.str.upper() + ms.reference_allele.str.upper()
    ms["palindromic"] = pal.isin(PALINDROMIC)
    n_pal = int(ms.palindromic.sum())
    ms_ok = ms[~ms.palindromic].copy()
    log(f"    palindromic (A/T, C/G) removed: {n_pal:,} "
        f"({100 * n_pal / len(ms):.2f}% of mapped) -> {len(ms_ok):,} instruments")

    ms_ok.rename(columns={"chr": "chr_hg19"})[
        ["cpg", "snp", "chr_hg19", "snp_pos", "pos38", "reference_allele",
         "other_allele", "eaf", "beta", "se", "z", "neg_log10p", "n_studies"]
    ].to_csv(os.path.join(a.out, "msites_hg38.tsv"), sep="\t", index=False)

    target = ms_ok[["cpg", "chrom", "pos38", "reference_allele", "other_allele",
                    "eaf", "beta", "se", "z", "p", "neg_log10p"]].copy()
    target["pos38"] = target.pos38.astype(np.int64)
    n_instruments = len(target)
    log(f"    instruments entering the screen: {n_instruments:,}")

    stats["pace"] = {
        "msites_total": int(n0),
        "autosomal_instruments": int(n_attempt),
        "lifted_hg38": int(len(ms)),
        "lift_unmapped": int(n_unmapped),
        "palindromic_excluded": int(n_pal),
        "instruments_tested": int(n_instruments),
    }

    # ------------------------------------------------------- screen outcomes
    for name, fname, kind, cases, controls in OUTCOMES:
        log(f"== 2. {name} ==")
        g, gs = read_gwas(os.path.join(a.raw, fname), kind)
        m = g.merge(target, on=["chrom", "pos38"], how="inner")
        n_match = len(m)
        log(f"    variants read {gs['variants_read']:,}; autosomal "
            f"{gs['autosomal_variants']:,}; instrument-position matches {n_match:,}")
        if m.empty:
            log("    no matches - check that the GWAS file uses GRCh38 positions")
            return 1

        # orient GWAS effects to the PACE effect allele
        fwd = (m.gwas_eff == m.reference_allele) & (m.gwas_oth == m.other_allele)
        rev = (m.gwas_eff == m.other_allele) & (m.gwas_oth == m.reference_allele)
        m = m[fwd | rev].copy()
        m["gwas_beta"] = np.where(m.gwas_eff == m.reference_allele,
                                  m.gwas_beta, -m.gwas_beta)
        m["gwas_af"] = np.where(m.gwas_eff == m.reference_allele,
                                m.af_eff, 1.0 - m.af_eff)
        n_allele = len(m)
        r = float(np.corrcoef(m.eaf, m.gwas_af)[0, 1])
        log(f"    allele-consistent {n_allele:,}; PACE eaf vs GWAS af r = {r:.3f}")

        mp = m.other_allele.str.upper() + m.reference_allele.str.upper()
        m["palindromic"] = mp.isin(PALINDROMIC)
        n_pal_g = int(m.palindromic.sum())
        m = m[~m.palindromic].copy()
        n_dup = int(m.duplicated("cpg").sum())
        m = m.drop_duplicates("cpg").copy()
        log(f"    palindromic after matching {n_pal_g:,}; duplicate CpG rows "
            f"{n_dup:,}; unique CpG {len(m):,}")

        d = smr_estimate(
            m.rename(columns={"p": "pace_p"}), n_instruments)
        d = d.rename(columns={"reference_allele": "pace_effect_allele",
                              "other_allele": "pace_non_effect_allele",
                              "eaf": "eaf_effect"})
        d.to_csv(os.path.join(a.out, f"smr_{name}.tsv"), sep="\t",
                 index=False, float_format="%.10g")
        hits = d[d.smr_hit]
        hits.to_csv(os.path.join(a.out, f"smr_hits_{name}.tsv"), sep="\t",
                    index=False, float_format="%.10g")

        thr = float(d.bonferroni_thr.iloc[0])
        log(f"    tested {len(d):,}; Bonferroni {thr:.4e}; SMR hits {len(hits):,}")
        for _, h in hits.iterrows():
            log(f"      {h.cpg}  {h.chrom}:{h.pos38}  b_SMR={h.b_smr:+.4f} "
                f"(se {h.se_smr:.4f})  P_SMR={h.p_smr:.2e}")

        stats[name] = {**gs, "position_matches": int(n_match),
                       "allele_consistent": int(n_allele),
                       "palindromic_dropped": n_pal_g,
                       "duplicate_cpg_removed": int(n_dup),
                       "eaf_af_r": round(r, 4), "cases": cases,
                       "controls": controls, "tested": int(len(d)),
                       "bonferroni_thr": thr, "smr_hits": int(len(hits)),
                       "hits": hits[["cpg", "chrom", "pos38", "b_smr",
                                     "se_smr", "p_smr"]].to_dict("records")}

    pe_hits = {h["cpg"] for h in stats["finngen_pe"]["hits"]}
    nc_hits = {h["cpg"] for h in stats["finngen_negctrl"]["hits"]}
    rc_hits = {h["cpg"] for h in stats["gcst90301704"]["hits"]}
    stats["cross"] = {
        "pe_hits": sorted(pe_hits),
        "pe_and_negctrl": sorted(pe_hits & nc_hits),
        "pe_and_gcst": sorted(pe_hits & rc_hits),
        "pe_only": sorted(pe_hits - nc_hits),
    }
    with open(os.path.join(a.out, "stats.json"), "w") as f:
        json.dump(stats, f, indent=2)

    # --------------------------------------------------------- expectations
    got = {
        "msites_total": stats["pace"]["msites_total"],
        "lifted_hg38": stats["pace"]["lifted_hg38"],
        "lift_unmapped": stats["pace"]["lift_unmapped"],
        "instruments_tested": n_instruments,
        "tested:finngen_pe": stats["finngen_pe"]["tested"],
    }
    log("\n== 3. counts versus the paper ==")
    bad = 0
    for k, want in EXPECTED.items():
        have = got[k]
        flag = "ok" if have == want else "DIFF"
        if have != want:
            bad += 1
        log(f"    {k:<20} expected {want:>9,}   obtained {have:>9,}   {flag}")

    if bad and a.strict:
        log("\nFAIL (--strict): at least one headline count differs from the paper.\n"
            "Intermediate counts depend on the UCSC chain release and on the exact\n"
            "PACE/GWAS file versions; see data/README.md.")
        return 1
    if bad:
        log("\nNote: intermediate counts depend on the chain release and on the\n"
            "PACE/GWAS file versions in data/raw. Re-run with --strict to make this\n"
            "fatal.")
    log("== done ==")
    return 0


if __name__ == "__main__":
    sys.exit(main())
