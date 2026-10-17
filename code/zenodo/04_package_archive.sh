#!/usr/bin/env bash
# 04_package_archive.sh -- zip data/ for the Zenodo record
# ------------------------------------------------------------------------------
# Run from the repository root after the pipeline and 01-03_*.R:
#   bash code/zenodo/04_package_archive.sh [output_dir]
# Writes four archives plus a manifest (path, bytes, md5) for each:
#   tree-methanogens_raw.zip        data/raw/ except picrust/
#   tree-methanogens_picrust.zip    data/raw/picrust/ (PICRUSt2 outputs, ~7 GB unzipped; read only by
#                                   stage A 01_import/11_picrust_pathways.R, whose small product is in
#                                   the processed archive)
#   tree-methanogens_processed.zip  data/processed/
#   tree-methanogens_compiled.zip   data/compiled/
# data/deprecated/ and .DS_Store files are excluded. Default output directory is
# ~/zenodo_upload/tree-methanogens (outside the repository and any synced folder).
# ------------------------------------------------------------------------------
set -euo pipefail
out="${1:-$HOME/zenodo_upload/tree-methanogens}"
mkdir -p "$out"
[ -d data/raw ] || { echo "run from the repository root" >&2; exit 1; }

manifest() {  # manifest <zip name> <find args...>
  local name="$1"; shift
  echo "path,bytes,md5" > "$out/MANIFEST_${name}.csv"
  find "$@" -type f ! -name .DS_Store ! -name .gitkeep -print0 | sort -z |
    while IFS= read -r -d '' f; do
      printf '"%s",%s,%s\n' "$f" "$(stat -f %z "$f")" "$(md5 -q "$f")"
    done >> "$out/MANIFEST_${name}.csv"
}

rm -f "$out"/tree-methanogens_*.zip
zip -qr "$out/tree-methanogens_raw.zip" data/raw -x "data/raw/picrust/*" "*.DS_Store" "*/.gitkeep"
manifest raw data/raw -path data/raw/picrust -prune -o
zip -qr "$out/tree-methanogens_picrust.zip" data/raw/picrust -x "*.DS_Store" "*/.gitkeep"
manifest picrust data/raw/picrust
zip -qr "$out/tree-methanogens_processed.zip" data/processed -x "*.DS_Store" "*/.gitkeep"
manifest processed data/processed
zip -qr "$out/tree-methanogens_compiled.zip" data/compiled -x "*.DS_Store" "*/.gitkeep"
manifest compiled data/compiled
ls -lh "$out"
