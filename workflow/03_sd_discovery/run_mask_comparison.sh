#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
project_dir="$(cd "$analysis_dir/.." && pwd)"
biser_bin="${BISER_BIN:-$(command -v biser || true)}"
source_fa="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.source_softmasked.fa"
manual_fa="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
result_root="$analysis_dir/runs"
statistics_dir="$analysis_dir/statistics"
threads=6

if [[ -z "$biser_bin" ]]; then
  printf 'ERROR: biser was not found. Add it to PATH or set BISER_BIN.\n' >&2
  exit 127
fi

mkdir -p "$result_root/source_softmask" "$result_root/annotation_extended_softmask" "$statistics_dir"

"$biser_bin" --version > "$statistics_dir/tool_version.txt"
printf 'threads\t%d\nparameters\tdefault thresholds; -t %d -o OUTPUT FASTA\n' \
  "$threads" "$threads" > "$statistics_dir/run_parameters.tsv"

run_one() {
  local label="$1"
  local fasta="$2"
  local out_dir="$result_root/$label"

  printf '%q -t %q -o %q %q\n' \
    "$biser_bin" "$threads" "$out_dir/biser_out" "$fasta" \
    > "$out_dir/command.txt"

  /usr/bin/time -v "$biser_bin" \
    -t "$threads" \
    -o "$out_dir/biser_out" \
    "$fasta" \
    > "$out_dir/stdout.log" \
    2> "$out_dir/resource_usage.log"
}

# Run serially so both conditions receive the same CPU allocation and do not
# compete for memory or I/O.
run_one source_softmask "$source_fa"
run_one annotation_extended_softmask "$manual_fa"
