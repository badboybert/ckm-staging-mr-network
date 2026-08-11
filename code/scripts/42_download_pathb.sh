#!/bin/bash
# --- portable roots (injected by scripts/build_repo.py) ---
P4_BASE="${CKM_P4_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
P4_ROOT="${CKM_P4_ROOT:-$(dirname "$P4_BASE")}"
CKM_ROOT="${CKM_ROOT:-$(dirname "$P4_ROOT")}"
# ----------------------------------------------------------
# Step 42: download the remaining Path-B GWAS (HF subtypes + WHRadjBMI), verifying md5 per file.
set -u
RAW="$P4_BASE/data/raw"
cd "$RAW" || exit 1
GWAS="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90728001-GCST90729000"

dl () {  # url  outname  expected_md5
  local url="$1" out="$2" md5="$3"
  echo "[dl] $out"
  curl -s -C - --retry 4 --retry-delay 5 -o "$out" "$url"
  local got; got=$(md5sum "$out" | cut -d' ' -f1)
  if [ "$got" = "$md5" ]; then echo "  OK  $(stat -c%s "$out") bytes  md5 $got"; else echo "  FAIL md5 $got != $md5"; fi
}

dl "$GWAS/GCST90728696/GCST90728696.tsv.gz" "HF_niHF_GCST90728696_EUR.tsv.gz"     "808b36f8284a6b7f80af809632b703df"
dl "$GWAS/GCST90728697/GCST90728697.tsv.gz" "HF_niHFpEF_GCST90728697_EUR.tsv.gz"  "c3bb00d4e4c462cb2d2c01f0f3ec3d84"
dl "$GWAS/GCST90728698/GCST90728698.tsv.gz" "HF_niHFrEF_GCST90728698_EUR.tsv.gz"  "ac57cfc846354e41af6600e50104b0e8"
dl "https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz" \
   "WHRadjBMI_Pulit2019_EUR.txt.gz" "500ee10b8406dc80394ad8cc968a4eb1"
echo "[done] Path-B downloads complete"
