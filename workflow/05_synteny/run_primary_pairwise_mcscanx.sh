#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
analysis_dir=$(cd "$script_dir/.." && pwd)
inputs_dir="$analysis_dir/inputs"
pair_inputs_dir="$analysis_dir/pair_inputs"
results_dir="$analysis_dir/results"
logs_dir="$analysis_dir/logs"
db_dir="$results_dir/reference_blast_db"

threads=${THREADS:-8}
mcscanx="${MCSCANX_BIN:-$(command -v MCScanX || true)}"
blastp_bin="${BLASTP_BIN:-$(command -v blastp || true)}"
makeblastdb_bin="${MAKEBLASTDB_BIN:-$(command -v makeblastdb || true)}"
codes=(Alyrata Bstricta Dstrictus Cviolacea)

if [[ -z "$mcscanx" || -z "$blastp_bin" || -z "$makeblastdb_bin" ]]; then
  printf 'ERROR: MCScanX, blastp, and makeblastdb are required; set the corresponding *_BIN variables if needed.\n' >&2
  exit 127
fi

mkdir -p "$results_dir" "$logs_dir" "$db_dir"

{
  printf 'field\tvalue\n'
  printf 'run_started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'mcscanx\t%s\n' "$mcscanx"
  printf 'mcscanx_sha256\t%s\n' "$(sha256sum "$mcscanx" | cut -d' ' -f1)"
  printf 'blastp\t%s\n' "$blastp_bin"
  printf 'makeblastdb\t%s\n' "$makeblastdb_bin"
  printf 'threads\t%s\n' "$threads"
  printf 'blast_evalue\t1e-10\n'
  printf 'blast_max_alignments\t5\n'
  printf 'mcscanx_match_size_s\t3\n'
  printf 'mcscanx_max_gaps_m\t2\n'
  printf 'mcscanx_overlap_window_w\t0\n'
} > "$logs_dir/run_environment.primary_jgi_20260724.tsv"

"$makeblastdb_bin" -in "$inputs_dir/Atha.pep" -dbtype prot -parse_seqids \
  -out "$db_dir/Atha" > "$logs_dir/makeblastdb.log" 2>&1

printf 'pair\tblast_command\tmcscanx_command\n' \
  > "$logs_dir/commands.primary_jgi_20260724.tsv"
for code in "${codes[@]}"; do
  pair="Atha_${code}"
  pair_dir="$results_dir/$pair"
  prefix="$pair_dir/$pair"
  mkdir -p "$pair_dir"
  cp "$pair_inputs_dir/$pair.gff" "$prefix.gff"

  blast_command="$blastp_bin -query $inputs_dir/$code.pep -db $db_dir/Atha -out $prefix.blast -evalue 1e-10 -num_threads $threads -outfmt 6 -num_alignments 5"
  mcscanx_command="$mcscanx -s 3 -m 2 -w 0 $pair"
  printf '%s\t%s\t%s\n' "$pair" "$blast_command" "$mcscanx_command" \
    >> "$logs_dir/commands.primary_jgi_20260724.tsv"

  echo "[$(date --iso-8601=seconds)] Running $pair"
  /usr/bin/time -v "$blastp_bin" \
    -query "$inputs_dir/$code.pep" \
    -db "$db_dir/Atha" \
    -out "$prefix.blast" \
    -evalue 1e-10 \
    -num_threads "$threads" \
    -outfmt 6 \
    -num_alignments 5 \
    > "$logs_dir/$pair.blast.stdout.log" \
    2> "$logs_dir/$pair.blast.time.log"
  (
    cd "$pair_dir"
    /usr/bin/time -v "$mcscanx" -s 3 -m 2 -w 0 "$pair"
  ) > "$logs_dir/$pair.mcscanx.stdout.log" \
    2> "$logs_dir/$pair.mcscanx.time.log"
done

printf 'run_finished\t%s\n' "$(date --iso-8601=seconds)" \
  >> "$logs_dir/run_environment.primary_jgi_20260724.tsv"
