#!/usr/bin/env bash
set -euo pipefail

project_dir="${SD_PROJECT_ROOT:?请设置 SD_PROJECT_ROOT 指向分析数据目录}"
biser_bin="${BISER_BIN:-$(command -v biser || true)}"
samtools_bin="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
threads="${BISER_THREADS:-6}"
rm_input="${REPEATMASKER_FASTA:-$project_dir/01_reference/repeatmasker_arabidopsis/GCA_978657495.1_TAIR12_genomic.uppercase.fna.masked}"
prepared_fa="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.repeatmasker_arabidopsis_softmasked.fa"
out_dir="$project_dir/03_biser_segmental_duplication/runs/repeatmasker_arabidopsis_softmask"

if [[ -z "$biser_bin" ]]; then
  printf 'ERROR: biser 不在 PATH 中，请设置 BISER_BIN。\n' >&2
  exit 127
fi
if [[ -z "$samtools_bin" ]]; then
  printf 'ERROR: samtools 不在 PATH 中，请设置 SAMTOOLS_BIN。\n' >&2
  exit 127
fi
if [[ ! -s "$rm_input" ]]; then
  printf 'ERROR: RepeatMasker soft-mask FASTA 不存在或为空：%s\n' "$rm_input" >&2
  exit 1
fi
if [[ -e "$out_dir/biser_out" || -e "$out_dir/biser_out.elem.txt" ]]; then
  printf 'ERROR: %s 已有 BISER 结果，请更换输出目录或先保留现有结果。\n' "$out_dir" >&2
  exit 1
fi

mkdir -p "$(dirname "$prepared_fa")" "$out_dir"

awk '
  BEGIN {
    map["OZ408683.1"]="Chr1"; map["OZ408684.1"]="Chr2";
    map["OZ408685.1"]="Chr3"; map["OZ408686.1"]="Chr4";
    map["OZ408687.1"]="Chr5"
  }
  /^>/ {
    name=substr($1,2)
    if (!(name in map)) {
      print "ERROR: unexpected FASTA record: " name > "/dev/stderr"
      exit 2
    }
    print ">" map[name]
    next
  }
  {print}
' "$rm_input" > "$prepared_fa"

"$samtools_bin" faidx "$prepared_fa"

read -r total_bp lower_bp n_bp records < <(
  awk '
    /^>/ {records++; next}
    {
      total += length($0)
      for (i=1; i<=length($0); i++) {
        b=substr($0,i,1)
        if (b ~ /[acgt]/) lower++
        else if (b ~ /[Nn]/) n++
      }
    }
    END {print total+0, lower+0, n+0, records+0}
  ' "$prepared_fa"
)

if [[ "$total_bp" -ne 142481245 || "$lower_bp" -ne 37623317 || "$n_bp" -ne 0 || "$records" -ne 5 ]]; then
  printf 'ERROR: FASTA QC 不符：total=%s lower=%s N=%s records=%s\n' \
    "$total_bp" "$lower_bp" "$n_bp" "$records" >&2
  exit 1
fi

"$biser_bin" --version > "$out_dir/tool_version.txt"
printf 'condition\trepeatmasker_arabidopsis_softmask\nthreads\t%s\nparameters\tdefault thresholds; -t %s -o OUTPUT FASTA\ninput_total_bp\t%s\ninput_lowercase_bp\t%s\ninput_N_bp\t%s\n' \
  "$threads" "$threads" "$total_bp" "$lower_bp" "$n_bp" \
  > "$out_dir/run_parameters.tsv"
printf '%q -t %q -o %q %q\n' \
  "$biser_bin" "$threads" "$out_dir/biser_out" "$prepared_fa" \
  > "$out_dir/command.txt"

/usr/bin/time -v "$biser_bin" \
  -t "$threads" -o "$out_dir/biser_out" "$prepared_fa" \
  > "$out_dir/stdout.log" 2> "$out_dir/resource_usage.log"

printf 'BISER 完成：%s\n' "$out_dir/biser_out"
