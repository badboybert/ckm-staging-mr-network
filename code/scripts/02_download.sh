#!/usr/bin/env bash
# --- portable roots (injected by scripts/build_repo.py) ---
P4_BASE="${CKM_P4_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
P4_ROOT="${CKM_P4_ROOT:-$(dirname "$P4_BASE")}"
CKM_ROOT="${CKM_ROOT:-$(dirname "$P4_ROOT")}"
# ----------------------------------------------------------
# Step 2: resumable, integrity-checked download of missing EUR exposure GWAS.
# Skips files already complete (remote Content-Length == local size). Logs progress.
set -u
RAW="$P4_BASE/data/raw"
LOG="$P4_BASE/logs/download.log"
mkdir -p "$RAW"
echo "=== download start $(date) ===" >> "$LOG"

# name|url|outfile
GC="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
GLGC="https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific"
TARGETS=(
  "BMI|$GC/GCST006001-GCST007000/GCST006900/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz|BMI_Yengo2018_EUR.txt.gz"
  "SBP|$GC/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz|SBP_Evangelou2018_EUR.txt.gz"
  "HbA1c|$GC/GCST90014001-GCST90015000/GCST90014006/harmonised/34017140-GCST90014006-EFO_0004541.h.tsv.gz|HbA1c_MAGIC_EUR.h.tsv.gz"
  "FI|$GC/GCST90002001-GCST90003000/GCST90002238/harmonised/34059833-GCST90002238-EFO_0004467-Build37.f.tsv.gz|FI_MAGIC_EUR.f.tsv.gz"
  "HDL|$GLGC/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz|HDL_GLGC2021_EUR.gz"
  "TC|$GLGC/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz|TC_GLGC2021_EUR.gz"
  "TG|$GLGC/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz|TG_GLGC2021_EUR.gz"
)
for t in "${TARGETS[@]}"; do
  IFS='|' read -r name url out <<< "$t"
  dest="$RAW/$out"
  remote=$(curl -sI -m 30 "$url" 2>/dev/null | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
  local=0; [ -f "$dest" ] && local=$(stat -c%s "$dest" 2>/dev/null || echo 0)
  if [ -n "$remote" ] && [ "$local" = "$remote" ]; then
    echo "[$(date +%H:%M:%S)] SKIP $name (complete, $local B)" | tee -a "$LOG"; continue
  fi
  echo "[$(date +%H:%M:%S)] GET  $name  remote=${remote:-?}B local=${local}B -> $out" | tee -a "$LOG"
  curl -s -m 3600 -C - -o "$dest" "$url" 2>>"$LOG"
  new=$(stat -c%s "$dest" 2>/dev/null || echo 0)
  ok="?"; if [ -n "$remote" ]; then [ "$new" = "$remote" ] && ok="OK" || ok="SIZE_MISMATCH"; fi
  gz="?"; gzip -t "$dest" 2>/dev/null && gz="gz_ok" || gz="gz_BAD"
  echo "[$(date +%H:%M:%S)] DONE $name  size=$new  $ok  $gz" | tee -a "$LOG"
done
echo "=== download end $(date) ===" >> "$LOG"
