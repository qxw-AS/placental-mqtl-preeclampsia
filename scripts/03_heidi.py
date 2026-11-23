#!/usr/bin/env python3
"""Step 3 - build BESD inputs and run HEIDI / SMR-multi for the SMR hits.

Step 2 screens all CpGs with the ratio estimator in Python. This step hands the
small set of significant CpGs to the SMR binary, which needs a BESD store and
variant identifiers, and is the only stage where rsIDs are required.

  work/smr_hits_finngen_pe.tsv     CpGs carried forward from step 2
  data/raw/pace_cis_mqtl_top_snps.tsv
  data/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr>.bim   rsID lookup by position
  data/ref/cpg_annotation.tsv      cpg, chr, pos_hg19 (EPIC manifest coordinates)
  data/raw/finngen_R13_*.gz, GCST90301704.h.tsv.gz

  -> work/smr_input/<cpg>.esd, pace_hits.flist
  -> work/smr_input/<outcome>.ma
  -> work/smr_commands.sh          the SMR v1.3.1 invocations
  -> work/smr_input/smr_multi_<outcome>.msmr   (only with --run)

The CpG coordinate does not enter the SMR ratio itself; it defines the probe
position written to the BESD store and therefore the cis window used by
SMR-multi.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

GWAS = [
    ("finngen_pe", "finngen_R13_O15_PREECLAMPS.gz", "finngen", 9151, 258775),
    ("finngen_negctrl", "finngen_R13_O15_PREEC_OR_FETGRO.gz", "finngen", 14616, 267293),
    ("gcst90301704", "GCST90301704.h.tsv.gz", "gcst", 1728, 192399),
]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def load_bim(bim_dir: str, chroms) -> dict:
    """(pos, frozenset(alleles)) -> rsID, per chromosome."""
    maps = {}
    for c in sorted(chroms):
        path = os.path.join(bim_dir, f"1000G.EUR.QC.{c}.bim")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path}\n"
                "The 1000G EUR plink .bim files are shipped with the SMR reference\n"
                "data (1000G_EUR_Phase3_plink). Put them under data/ref/.")
        m = {}
        with open(path) as fh:
            for line in fh:
                a = line.split()
                # keep the first rsID seen at a position/allele combination
                m.setdefault((int(a[3]), frozenset((a[4], a[5]))), a[1])
        maps[str(c)] = m
        log(f"    bim chr{c}: {len(m):,} positions")
    return maps


def write_ma(path: str, out: str, kind: str, rsids: set, a1: pd.DataFrame,
             cases: int, controls: int) -> int:
    """Write an SMR .ma file restricted to the instrument rsIDs."""
    if kind == "finngen":
        cols = {"#chrom": str, "rsids": str, "ref": str, "alt": str,
                "beta": np.float64, "sebeta": np.float64,
                "pval": np.float64, "af_alt": np.float64}
        reader = pd.read_csv(path, sep="\t", usecols=list(cols), dtype=cols,
                             chunksize=2_000_000, low_memory=False)
    else:
        cols = {"rsid": str, "effect_allele": str, "other_allele": str,
                "effect_allele_frequency": np.float64, "beta": np.float64,
                "standard_error": np.float64, "p_value": np.float64}
        reader = pd.read_csv(path, sep="\t", usecols=list(cols), dtype=cols,
                             chunksize=2_000_000, low_memory=False)

    keep, n = [], 0
    for chunk in reader:
        n += len(chunk)
        if kind == "finngen":
            ex = chunk[chunk.rsids.notna()].copy()
            if ex.empty:
                continue
            ex["rs"] = ex.rsids.str.split(",")
            ex = ex.explode("rs")
            ex = ex[ex.rs.isin(rsids)]
        else:
            ex = chunk[chunk.rsid.isin(rsids)].rename(columns={"rsid": "rs"})
        keep.append(ex)

    if not keep:
        raise RuntimeError(f"no instrument rsIDs found in {path}")
    k = pd.concat(keep, ignore_index=True)
    k = k.drop_duplicates("rs", keep="first")
    if k.empty:
        raise RuntimeError(f"no instrument rsIDs found in {path}")

    if kind == "finngen":
        eff, oth, freq, beta, se, p = (k.alt, k.ref, k.af_alt, k.beta,
                                       k.sebeta, k.pval)
    else:
        eff, oth, freq, beta, se, p = (k.effect_allele, k.other_allele,
                                       k.effect_allele_frequency, k.beta,
                                       k.standard_error, k.p_value)

    ref = a1.loc[k.rs, "reference_allele"].to_numpy()
    alt = a1.loc[k.rs, "other_allele"].to_numpy()
    fwd = (eff.to_numpy() == ref) & (oth.to_numpy() == alt)
    rev = (eff.to_numpy() == alt) & (oth.to_numpy() == ref)
    ok = fwd | rev
    if not ok.all():
        log(f"    dropped {int((~ok).sum())} rsIDs with non-matching alleles")

    o = pd.DataFrame({
        "SNP": k.rs.to_numpy()[ok],
        "A1": ref[ok],
        "A2": alt[ok],
        "freq": np.where(fwd[ok], freq.to_numpy()[ok], 1 - freq.to_numpy()[ok]),
        "b": np.where(fwd[ok], beta.to_numpy()[ok], -beta.to_numpy()[ok]),
        "se": se.to_numpy()[ok],
        "p": p.to_numpy()[ok],
        "n": cases + controls,
    })
    o.to_csv(out, sep="\t", index=False, float_format="%.6g")
    log(f"    variants {n:,}; matched rsIDs {len(k):,}/{len(rsids):,}; "
        f"allele-consistent {len(o):,} -> {out}")
    return len(o)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", default="work/smr_hits_finngen_pe.tsv")
    ap.add_argument("--pace", default="data/raw/pace_cis_mqtl_top_snps.tsv")
    ap.add_argument("--raw", default="data/raw")
    ap.add_argument("--bim-dir", default="data/ref/1000G_EUR_Phase3_plink")
    ap.add_argument("--cpg-pos", default="data/ref/cpg_annotation.tsv")
    ap.add_argument("--out", default="work/smr_input")
    ap.add_argument("--smr-bin", default="bin/smr")
    ap.add_argument("--bfile", default="data/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC",
                    help="plink prefix or list file passed to SMR --bfile for LD")
    ap.add_argument("--run", action="store_true",
                    help="execute the SMR commands after writing the inputs")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)

    hits = pd.read_csv(a.hits, sep="\t")
    HITS = hits.cpg.tolist()
    log(f"== 1. CpGs carried forward from step 2: {len(HITS)} ==")
    for c in HITS:
        log(f"    {c}")

    cpg = pd.read_csv(a.cpg_pos, sep="\t")
    cpg["cpg"] = cpg.cpg.astype(str)
    missing = set(HITS) - set(cpg.cpg)
    if missing:
        log(f"    missing CpG coordinates in {a.cpg_pos}: {sorted(missing)}")
        return 1
    cpg_pos = cpg[cpg.cpg.isin(HITS)].drop_duplicates("cpg").set_index("cpg")[
        ["chr", "pos_hg19"]].to_dict("index")

    # --------------------------------------------------- PACE pairs -> rsIDs
    log("== 2. map PACE SNPs to 1000G rsIDs ==")
    pace = pd.read_csv(a.pace, sep="\t", dtype=str, low_memory=False)
    d = pace[pace.cpg.isin(HITS)].copy()
    sp = d.snp.str.split(":", expand=True)
    d["chr"] = sp[0].astype(int)
    d["pos19"] = sp[1].astype(int)
    for c in ["eaf", "beta", "se", "p-value"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["eaf", "beta", "se"]).copy()

    maps = load_bim(a.bim_dir, d.chr.unique())
    d["rsid"] = [maps[str(c)].get((pos, frozenset((a1, a2))))
                 for c, pos, a1, a2 in zip(d.chr, d.pos19,
                                           d.reference_allele, d.other_allele)]
    log(f"    pairs {len(d):,}; mapped {d.rsid.notna().sum():,}; "
        f"unique rsIDs {d.rsid.nunique():,}")
    d = d[d.rsid.notna()].copy()

    eff = d.groupby("rsid")[["reference_allele", "other_allele"]].nunique()
    assert (eff <= 1).all().all(), eff[(eff > 1).any(axis=1)]
    a1 = d.drop_duplicates("rsid").set_index("rsid")[
        ["reference_allele", "other_allele"]]
    log(f"    unique instruments: {len(a1):,}")

    # ----------------------------------------------------------- ESD + flist
    log("== 3. write ESD + flist ==")
    flist = ["Chr ProbeID GeneticDistance ProbeBp Gene Orientation PathOfEsd"]
    for c in HITS:
        s = d[d.cpg == c].sort_values("pos19")
        esd = s[["chr", "rsid", "pos19", "reference_allele", "other_allele",
                 "eaf", "beta", "se", "p-value"]].copy()
        esd.columns = ["Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p"]
        path = os.path.abspath(os.path.join(a.out, f"{c}.esd"))
        esd.to_csv(path, sep="\t", index=False, float_format="%.6g")
        cc, bb = cpg_pos[c]["chr"], cpg_pos[c]["pos_hg19"]
        flist.append(f"{int(cc)}\t{c}\t0\t{int(bb)}\tNA\tNA\t{path}")
        log(f"    {c}: {len(s):,} SNPs -> {path}")
    flist_path = os.path.join(a.out, "pace_hits.flist")
    with open(flist_path, "w") as fh:
        fh.write("\n".join(flist) + "\n")

    # ------------------------------------------------------------- .ma files
    log("== 4. GWAS .ma files ==")
    RSSET = set(a1.index)
    for name, fname, kind, cases, controls in GWAS:
        write_ma(os.path.join(a.raw, fname), os.path.join(a.out, f"{name}.ma"),
                 kind, RSSET, a1, cases, controls)

    # ------------------------------------------------------------ SMR calls
    besd = os.path.join(a.out, "pace_hits")
    cmds = [f"{a.smr_bin} --make-besd --eqtl-flist {flist_path} --out {besd}"]
    for name, *_ in GWAS:
        cmds.append(
            f"{a.smr_bin} --beqtl-summary {besd}.besd "
            f"--gwas-summary {os.path.join(a.out, name + '.ma')} "
            f"--multi --bfile {a.bfile} "
            f"--out {os.path.join(a.out, 'smr_multi_' + name)}")
    script = os.path.join("work", "smr_commands.sh")
    os.makedirs("work", exist_ok=True)
    with open(script, "w", newline="\n") as fh:
        fh.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        fh.write("\n".join(cmds) + "\n")
    log(f"== 5. SMR v1.3.1 invocations written to {script} ==")
    for c in cmds:
        log("    " + c)

    if not a.run:
        log("\nInputs are complete. Re-run with --run to execute the SMR binary.")
        return 0

    if not (os.path.exists(a.smr_bin) or os.path.exists(a.smr_bin + ".exe")):
        log(f"    SMR binary not found at {a.smr_bin} (see environment/smr_version.txt)")
        return 2
    for c in cmds:
        log("$ " + c)
        subprocess.run(c, shell=True, check=True)
    log("== done ==")
    return 0


if __name__ == "__main__":
    sys.exit(main())
