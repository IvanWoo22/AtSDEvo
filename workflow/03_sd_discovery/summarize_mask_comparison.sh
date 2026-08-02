#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
project_dir="$(cd "$analysis_dir/.." && pwd)"
run_dir="$analysis_dir/runs"
statistics_dir="$analysis_dir/statistics"
added_mask="$project_dir/02_softmask_evaluation/intervals/newly_added_softmask.bed"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

printf 'condition\treported_sd_pairs\tdecomposition_elements\tputative_alignments\tintrachromosomal\tinterchromosomal\tmean_copy_length_bp\tmedian_copy_length_bp\tunion_copy_coverage_bp\tevents_overlapping_added_mask\tpct_events_overlapping_added_mask\twall_seconds\tmax_rss_kb\n' \
  > "$statistics_dir/summary.tsv"

for condition in source_softmask annotation_extended_softmask; do
  out="$run_dir/$condition/biser_out"
  elem="$run_dir/$condition/biser_out.elem.txt"
  stdout="$run_dir/$condition/stdout.log"
  time_log="$run_dir/$condition/resource_usage.log"

  pairs="$(wc -l < "$out")"
  elements="$(wc -l < "$elem")"
  putative="$(awk '/Total alignments:/{gsub(/,/,"",$3); print $3}' "$stdout")"
  read -r intra inter mean_len < <(
    awk -F '\t' '{a=$3-$2; b=$6-$5; sum+=a+b; if($1==$4)intra++; else inter++}
      END{printf "%d %d %.2f\n",intra,inter,sum/(2*NR)}' "$out"
  )
  awk -F '\t' '{print $3-$2; print $6-$5}' "$out" | LC_ALL=C sort -n \
    | awk '{v[NR]=$1} END{if(NR%2) print v[(NR+1)/2]; else print (v[NR/2]+v[NR/2+1])/2}' \
    > "$tmp_dir/median"
  median_len="$(cat "$tmp_dir/median")"

  awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$2,$3; print $4,$5,$6}' "$out" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - \
    | awk '{s+=$3-$2} END{print s+0}' > "$tmp_dir/coverage"
  union_coverage="$(cat "$tmp_dir/coverage")"

  awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,NR; print $4,$5,$6,NR}' "$out" \
    | bedtools intersect -a - -b "$added_mask" -u | cut -f4 | LC_ALL=C sort -u \
    | wc -l > "$tmp_dir/added_overlap_count"
  added_overlap_count="$(cat "$tmp_dir/added_overlap_count")"

  wall_raw="$(awk -F ': ' '/Elapsed \(wall clock\)/{print $2}' "$time_log")"
  wall_seconds="$(awk -v t="$wall_raw" 'BEGIN{n=split(t,a,":"); if(n==2)print a[1]*60+a[2]; else print a[1]*3600+a[2]*60+a[3]}')"
  max_rss="$(awk -F ': ' '/Maximum resident set size/{print $2}' "$time_log")"

  awk -v c="$condition" -v pairs="$pairs" -v elements="$elements" \
    -v putative="$putative" -v intra="$intra" -v inter="$inter" \
    -v mean="$mean_len" -v median="$median_len" -v coverage="$union_coverage" \
    -v overlap="$added_overlap_count" -v wall="$wall_seconds" -v rss="$max_rss" \
    'BEGIN{OFS="\t"; printf "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.1f\t%d\t%d\t%.4f\t%.2f\t%d\n",
      c,pairs,elements,putative,intra,inter,mean,median,coverage,overlap,100*overlap/pairs,wall,rss}' \
    >> "$statistics_dir/summary.tsv"

  cut -f1-6,9-10 "$out" | LC_ALL=C sort -u > "$tmp_dir/$condition.keys"

  awk -F '\t' -v condition="$condition" 'BEGIN{OFS="\t"}{a=$1; b=$4; if(a>b){t=a;a=b;b=t}; n[a"\t"b]++}
    END{for(k in n) print condition,k,n[k]}' "$out" \
    | LC_ALL=C sort -t $'\t' -k2,2 -k3,3 >> "$tmp_dir/chromosome_pairs"
done

{
  printf 'condition\tchromosome_1\tchromosome_2\tevents\n'
  cat "$tmp_dir/chromosome_pairs"
} > "$statistics_dir/chromosome_pair_counts.tsv"

shared="$(comm -12 "$tmp_dir/source_softmask.keys" "$tmp_dir/annotation_extended_softmask.keys" | wc -l)"
source_only="$(comm -23 "$tmp_dir/source_softmask.keys" "$tmp_dir/annotation_extended_softmask.keys" | wc -l)"
manual_only="$(comm -13 "$tmp_dir/source_softmask.keys" "$tmp_dir/annotation_extended_softmask.keys" | wc -l)"
printf 'comparison\tevents\nexact_coordinate_shared\t%s\nsource_only_exact\t%s\nannotation_extended_softmask_only_exact\t%s\n' \
  "$shared" "$source_only" "$manual_only" > "$statistics_dir/exact_coordinate_comparison.tsv"

awk -F '\t' 'NR==2{for(i=2;i<=NF;i++)source[i]=$i} NR==3{
    print "metric\tsource_original\tannotation_softmasked\tchange\tchange_pct";
    for(i=2;i<=NF;i++)
      printf "%s\t%s\t%s\t%.4f\t%.4f\n",header[i],source[i],$i,$i-source[i],100*($i-source[i])/source[i]
  } NR==1{for(i=2;i<=NF;i++)header[i]=$i}' "$statistics_dir/summary.tsv" \
  > "$statistics_dir/comparison_changes.tsv"
