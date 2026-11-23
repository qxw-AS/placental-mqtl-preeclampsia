#!/usr/bin/env bash
# Fetch the machine-retrievable public inputs into data/raw/ and data/ref/.
# Idempotent: existing non-empty files are skipped. Nothing here is committed to Git.
set -euo pipefail
cd "$(dirname "$0")/.."

RAW=data/raw
REF=data/ref
mkdir -p "$RAW" "$REF"

fetch () {  # fetch <url> <output-path>
  local url="$1" out="$2"
  if [ -s "$out" ]; then
    echo "[skip] $out already present"
  else
    echo "[get ] $out"
    curl -fL --retry 3 --retry-delay 5 -o "$out.part" "$url"
    mv "$out.part" "$out"
  fi
}

# --- 1. PACE placental cis-mQTL -------------------------------------------------
# The resource is a Streamlit app with no bulk endpoint and no versioned release.
# Export the top-SNP table by hand and save it here.
if [ ! -s "$RAW/pace_cis_mqtl_top_snps.tsv" ]; then
  cat <<'EOF'
[MANUAL] PACE top-SNP table must be exported from the web interface:
         https://pace-placenta-mqtl.streamlit.app   (version marker 2025-12-09)
         Save as data/raw/pace_cis_mqtl_top_snps.tsv (19 columns, tab-separated),
         then re-run this script.
EOF
fi

# --- 2. FinnGen R13 -------------------------------------------------------------
# Endpoint slugs and bucket path are those of the R13 public release.
FG_BASE="https://storage.googleapis.com/finngen-public-data-r13/summary_stats/release"
fetch "$FG_BASE/finngen_R13_O15_PREECLAMPS.gz"       "$RAW/finngen_R13_O15_PREECLAMPS.gz"
fetch "$FG_BASE/finngen_R13_O15_PREEC_OR_FETGRO.gz" "$RAW/finngen_R13_O15_PREEC_OR_FETGRO.gz"

# --- 3. GWAS Catalog GCST90301704 (harmonised, GRCh38) --------------------------
# The FTP tree groups accessions differently across releases, so try both layouts.
GC_OUT="$RAW/GCST90301704.h.tsv.gz"
if [ ! -s "$GC_OUT" ]; then
  for base in \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90301704" \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90301001-GCST90302000/GCST90301704"
  do
    if curl -fL --retry 2 -o "$GC_OUT.part" "$base/GCST90301704.h.tsv.gz"; then
      mv "$GC_OUT.part" "$GC_OUT"
      break
    fi
    rm -f "$GC_OUT.part"
  done
fi
if [ ! -s "$GC_OUT" ]; then
  echo "[MANUAL] GCST90301704 could not be fetched automatically."
  echo "         Download the harmonised file from https://www.ebi.ac.uk/gwas/studies/GCST90301704"
  echo "         and save it as $GC_OUT"
fi

# --- 4. UCSC liftOver chain -----------------------------------------------------
if [ ! -s "$REF/hg19ToHg38.over.chain" ] && [ ! -s "$REF/hg19ToHg38.over.chain.gz" ]; then
  fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
        "$REF/hg19ToHg38.over.chain.gz"
fi

# --- 5. 1000 Genomes Phase 3 EUR plink (.bim) ----------------------------------
# Ships with the SMR reference data (1000G_EUR_Phase3_plink); step 3 only needs the
# .bim files, and only for the chromosomes carrying hits.
if ! compgen -G "$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.*.bim" > /dev/null; then
  cat <<'EOF'
[MANUAL] Place the per-chromosome 1000G EUR .bim files at
         data/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr>.bim
         They are distributed with the SMR reference panel
         (https://cnsgenomics.com/software/smr/#DataResource).
EOF
fi

# --- 6. Record what was downloaded ---------------------------------------------
# File sizes and MD5s back the data-availability statement; keep them in sync with
# Table S5 of the manuscript.
echo "[check] recording size and MD5 of the downloaded inputs"
for f in "$RAW"/*.gz "$REF"/*.chain*; do
  [ -e "$f" ] || continue
  printf '%s  %s\n' "$(md5sum "$f" | cut -d' ' -f1)" "$f"
done

echo
echo "done. next: python scripts/01_msites.py"
