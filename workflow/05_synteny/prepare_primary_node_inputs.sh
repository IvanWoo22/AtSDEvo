#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
analysis_dir=$(cd "$script_dir/.." && pwd)
project_dir=$(cd "$analysis_dir/.." && pwd)
inputs_dir="$analysis_dir/inputs"
pair_inputs_dir="$analysis_dir/pair_inputs"
resources="$project_dir/04_outgroup_resources/primary_nodes_jgi_20260724/prepared"
python_bin="${PYTHON_BIN:-$(command -v python3 || true)}"
gffread_bin="${GFFREAD_BIN:-$(command -v gffread || true)}"

if [[ -z "$python_bin" || -z "$gffread_bin" ]]; then
  printf 'ERROR: python3 and gffread are required; set PYTHON_BIN/GFFREAD_BIN if needed.\n' >&2
  exit 127
fi

mkdir -p "$inputs_dir" "$pair_inputs_dir"

"$python_bin" \
  "$script_dir/extract_tair12_canonical_proteins.py" \
  --gff "$project_dir/01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3" \
  --genome "$project_dir/01_reference/prepared_data/TAIR12.Col-CC.source_softmasked.fa" \
  --outdir "$inputs_dir" \
  --gffread "$gffread_bin"

codes=(Alyrata Bstricta Dstrictus Cviolacea)
for code in "${codes[@]}"; do
  cp "$resources/$code.gff" "$inputs_dir/$code.gff"
  cp "$resources/$code.pep" "$inputs_dir/$code.pep"
  cp "$resources/$code.bed" "$inputs_dir/$code.bed"
  cp "$resources/$code.transcript_map.tsv" "$inputs_dir/$code.transcript_map.tsv"
  cp "$inputs_dir/Atha.gff" "$pair_inputs_dir/Atha_${code}.gff"
  sed -e '$a\' "$inputs_dir/$code.gff" >> "$pair_inputs_dir/Atha_${code}.gff"
done

"$python_bin" \
  "$script_dir/validate_mcscanx_inputs.py" --root "$analysis_dir"
