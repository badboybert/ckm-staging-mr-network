#!/usr/bin/env bash
# --- portable roots (injected by scripts/build_repo.py) ---
P4_BASE="${CKM_P4_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
P4_ROOT="${CKM_P4_ROOT:-$(dirname "$P4_BASE")}"
CKM_ROOT="${CKM_ROOT:-$(dirname "$P4_ROOT")}"
# ----------------------------------------------------------
set -u
RAW="$P4_BASE/data/raw"
LOG="$P4_BASE/logs/download_bbj.log"
BASE="https://humandbs.dbcls.jp/files/hum0014"
echo "=== BBJ QT download $(date) ===" >> "$LOG"
for tr in HbA1c HDL LDL TC TG; do
  dest="$RAW/${tr}_BBJ_hum0014.zip"
  if [ -f "$dest" ] && unzip -t "$dest" >/dev/null 2>&1; then echo "[skip] $tr ok" | tee -a "$LOG"; continue; fi
  echo "[get] $tr" | tee -a "$LOG"
  curl -s -m 1200 -C - -o "$dest" "$BASE/hum0014.v8.$tr.zip" 2>>"$LOG"
  if unzip -t "$dest" >/dev/null 2>&1; then
    unzip -o -q "$dest" -d "$RAW" && echo "[done] $tr $(stat -c%s "$dest")B unzipped" | tee -a "$LOG"
  else echo "[done] $tr ZIP_BAD" | tee -a "$LOG"; fi
done
echo "=== done $(date) ===" >> "$LOG"
