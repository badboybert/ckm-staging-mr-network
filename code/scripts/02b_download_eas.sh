#!/usr/bin/env bash
# --- portable roots (injected by scripts/build_repo.py) ---
P4_BASE="${CKM_P4_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
P4_ROOT="${CKM_P4_ROOT:-$(dirname "$P4_BASE")}"
CKM_ROOT="${CKM_ROOT:-$(dirname "$P4_ROOT")}"
# ----------------------------------------------------------
set -u
RAW="$P4_BASE/data/raw"
LOG="$P4_BASE/logs/download_eas.log"
GLGC="https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific"
echo "=== EAS lipid download $(date) ===" >> "$LOG"
declare -A T=(
 [HDL_GLGC2021_EAS.gz]="HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz"
 [LDL_GLGC2021_EAS.gz]="LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz"
 [TC_GLGC2021_EAS.gz]="TC_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz"
 [TG_GLGC2021_EAS.gz]="logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz"
)
for out in "${!T[@]}"; do
  url="$GLGC/${T[$out]}"; dest="$RAW/$out"
  remote=$(curl -sI -m 30 "$url" 2>/dev/null | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
  local=0; [ -f "$dest" ] && local=$(stat -c%s "$dest" 2>/dev/null||echo 0)
  if [ -n "$remote" ] && [ "$local" = "$remote" ]; then echo "[skip] $out done" | tee -a "$LOG"; continue; fi
  echo "[get] $out (${remote:-?}B)" | tee -a "$LOG"
  curl -s -m 3600 -C - -o "$dest" "$url" 2>>"$LOG"
  gzip -t "$dest" 2>/dev/null && echo "[done] $out $(stat -c%s "$dest")B gz_ok" | tee -a "$LOG" || echo "[done] $out gz_BAD" | tee -a "$LOG"
done
echo "=== done $(date) ===" >> "$LOG"
