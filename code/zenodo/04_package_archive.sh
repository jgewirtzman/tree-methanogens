#!/usr/bin/env bash
# 04_package_archive.sh -- build the Zenodo files (locally; uploading is done by hand)
# ------------------------------------------------------------------------------
# Run from the repository root after a full run (Rscript code/run_all.R --from raw):
#   bash code/zenodo/04_package_archive.sh [output_dir]
# Writes three archives, a manifest (path, bytes, md5) of each, and MD5SUMS.txt.
# Each data zip unpacks in place at the repository root.
#   tree-methanogens_data.zip     data/raw/ (without picrust/)   field, instrument and lab records
#                                 data/tables/                   the canonical tables
#                                 data/processed/                stage-A products (for --from B)
#                                 data/compiled/                 README_data.md and results/
#   tree-methanogens_picrust.zip  data/raw/picrust/ (PICRUSt2 predictions, ~7.9 GB unzipped).
#                                 Optional: read only by stage A's 01_import/06_picrust_pathways.R,
#                                 whose two small products are in the data zip.
#   tree-methanogens_code.zip     `git archive` of HEAD (code, renv.lock, README, folder layout),
#                                 unpacking into tree-methanogens/.
# Left out of the data zips:
#   data/deprecated/                                  retired files, not data of the paper
#   .DS_Store, Icon\r (Finder, Google Drive), .gitkeep (the code zip has them), ~$* (Office locks)
#   .RData, .Rhistory (R sessions; data/raw/inventory/spatial_data/ holds an unrelated one)
#   data/raw/ddpcr/Data Files/tree_ddpcr_full.R       a lab analysis script with personal paths
# The script stops if tracked files have uncommitted changes (so the code zip is what was
# run) or if data/compiled/ holds anything besides README_data.md and results/. Default
# output directory: ~/zenodo_upload/tree-methanogens (outside the repository).
# ------------------------------------------------------------------------------
set -euo pipefail
out="${1:-$HOME/zenodo_upload/tree-methanogens}"
[ -d data/raw ] && [ -f code/pipeline.csv ] || { echo "run from the repository root" >&2; exit 1; }
git diff --quiet HEAD -- || { echo "uncommitted changes to tracked files: commit first" >&2; exit 1; }
for d in data/raw/picrust data/processed data/tables data/compiled/results; do
  [ -d "$d" ] || { echo "missing $d -- run the pipeline first" >&2; exit 1; }
done
[ -f data/compiled/README_data.md ] || { echo "missing data/compiled/README_data.md" >&2; exit 1; }
extra="$(find data/compiled -mindepth 1 -maxdepth 1 ! -name README_data.md ! -name results ! -name .gitkeep)"
[ -z "$extra" ] || { echo "data/compiled/ holds more than README_data.md and results/:" >&2; echo "$extra" >&2; exit 1; }
mkdir -p "$out"
out="$(cd "$out" && pwd)"
commit="$(git rev-parse --short HEAD)"

# file size and md5, portable across GNU/Linux and macOS/BSD
size_of() { wc -c < "$1" | tr -d ' '; }
if command -v md5sum >/dev/null 2>&1; then md5_of() { md5sum "$1" | cut -d' ' -f1; }
else md5_of() { md5 -q "$1"; }; fi

# files that never go in the data archives (see the header)
JUNK=( -name .DS_Store -o -name .gitkeep -o -name $'Icon\r' -o -name '~$*'
       -o -name .RData -o -name .Rhistory
       -o -path 'data/raw/ddpcr/Data Files/tree_ddpcr_full.R' )
list_files() {  # list_files <dir>... : NUL-separated, sorted, junk and picrust/ removed
  find "$@" -path data/raw/picrust -prune -o -type f ! \( "${JUNK[@]}" \) -print0 | sort -z
}

pack() {  # pack <name> <dir>... : zip the listed files and write their manifest
  local name="$1"; shift
  local zip="$out/tree-methanogens_${name}.zip" man="$out/MANIFEST_${name}.csv"
  rm -f "$zip"
  echo "path,bytes,md5" > "$man"
  list_files "$@" | while IFS= read -r -d '' f; do
    printf '"%s",%s,%s\n' "$f" "$(size_of "$f")" "$(md5_of "$f")"
  done >> "$man"
  list_files "$@" | tr '\0' '\n' | zip -q -X "$zip" -@
  echo "  $name: $(($(wc -l < "$man") - 1)) files"
}

echo "Building the Zenodo files in $out (code at $commit)"
pack data data/raw data/tables data/processed data/compiled
# picrust/ is pruned by list_files; zip it on its own
rm -f "$out/tree-methanogens_picrust.zip"
echo "path,bytes,md5" > "$out/MANIFEST_picrust.csv"
find data/raw/picrust -type f ! \( "${JUNK[@]}" \) -print0 | sort -z | while IFS= read -r -d '' f; do
  printf '"%s",%s,%s\n' "$f" "$(size_of "$f")" "$(md5_of "$f")"
done >> "$out/MANIFEST_picrust.csv"
find data/raw/picrust -type f ! \( "${JUNK[@]}" \) | sort | zip -q -X "$out/tree-methanogens_picrust.zip" -@
echo "  picrust: $(($(wc -l < "$out/MANIFEST_picrust.csv") - 1)) files"
git archive --format=zip --prefix=tree-methanogens/ -o "$out/tree-methanogens_code.zip" HEAD
echo "  code: git archive of $commit"

( cd "$out" && for z in tree-methanogens_*.zip; do printf '%s  %s  %s bytes\n' "$(md5_of "$z")" "$z" "$(size_of "$z")"; done ) > "$out/MD5SUMS.txt"
echo "code commit: $(git rev-parse HEAD)" >> "$out/MD5SUMS.txt"
cat "$out/MD5SUMS.txt"
