#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
project_dir="$(cd "$analysis_dir/.." && pwd)"
source_fa="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.source_softmasked.fa"
final_fa="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
gff="$project_dir/01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3"
interval_dir="$analysis_dir/intervals"
statistics_dir="$analysis_dir/statistics"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

mkdir -p "$interval_dir" "$statistics_dir"

# Convert every maximal lowercase run in FASTA into a 0-based half-open BED
# interval without changing sequence names or loading the genome into memory.
fasta_lowercase_to_bed() {
  perl -ne '
    chomp;
    if (/^>(\S+)/) {
      $chr = $1;
      $offset = 0;
      next;
    }
    while (/([a-z]+)/g) {
      print join("\t", $chr, $offset + $-[1], $offset + $+[1]), "\n";
    }
    $offset += length($_);
  ' "$1"
}

fasta_lowercase_to_bed "$source_fa" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -i - > "$interval_dir/source_softmask.bed"

fasta_lowercase_to_bed "$final_fa" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -i - > "$interval_dir/final_softmask.bed"

# The newly masked set is exactly final minus source. -A is inappropriate here
# because partial overlaps must be trimmed rather than deleting whole records.
bedtools subtract \
  -a "$interval_dir/final_softmask.bed" \
  -b "$interval_dir/source_softmask.bed" \
  | bedtools merge -i - > "$interval_dir/newly_added_softmask.bed"

# One BED record per GFF feature. Type remains in column 4 for grouping.
awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ {print $1,$4-1,$5,$3}' "$gff" \
  | LC_ALL=C sort -k4,4 -k1,1 -k2,2n -k3,3n \
  > "$tmp_dir/features.by_type.bed"

bed_bp() {
  awk '{s += $3-$2} END{print s+0}' "$1"
}

overlap_bp() {
  bedtools intersect -a "$1" -b "$2" -wo \
    | awk '{s += $NF} END{print s+0}'
}

source_bp="$(bed_bp "$interval_dir/source_softmask.bed")"
final_bp="$(bed_bp "$interval_dir/final_softmask.bed")"
added_bp="$(bed_bp "$interval_dir/newly_added_softmask.bed")"

printf '%s\n' \
  $'feature_type\tfeature_records\tfeature_union_bp\tsource_sm_overlap_bp\tsource_sm_pct_of_feature\tsource_sm_pct_explained_by_type\tnew_sm_overlap_bp\tnew_sm_pct_of_feature\tnew_sm_pct_explained_by_type\tfinal_sm_overlap_bp\tfinal_sm_pct_of_feature\tselected_by_mask_script' \
  > "$statistics_dir/coverage_by_gff_type.tsv"

cut -f4 "$tmp_dir/features.by_type.bed" | LC_ALL=C sort -u \
  | while IFS= read -r feature_type; do
      awk -F '\t' -v t="$feature_type" 'BEGIN{OFS="\t"} $4==t{print $1,$2,$3}' \
        "$tmp_dir/features.by_type.bed" > "$tmp_dir/type.raw.bed"
      records="$(wc -l < "$tmp_dir/type.raw.bed")"
      LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$tmp_dir/type.raw.bed" \
        | bedtools merge -i - > "$tmp_dir/type.union.bed"
      feature_bp="$(bed_bp "$tmp_dir/type.union.bed")"
      source_overlap="$(overlap_bp "$tmp_dir/type.union.bed" "$interval_dir/source_softmask.bed")"
      added_overlap="$(overlap_bp "$tmp_dir/type.union.bed" "$interval_dir/newly_added_softmask.bed")"
      final_overlap="$(overlap_bp "$tmp_dir/type.union.bed" "$interval_dir/final_softmask.bed")"
      selected=no
      if [[ "$feature_type" == "mobile_genetic_element" || "${feature_type,,}" == *repeat* ]]; then
        selected=yes
      fi
      awk -v t="$feature_type" -v n="$records" -v f="$feature_bp" \
        -v so="$source_overlap" -v st="$source_bp" \
        -v ao="$added_overlap" -v at="$added_bp" \
        -v fo="$final_overlap" -v selected="$selected" 'BEGIN{OFS="\t";
          printf "%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%d\t%.4f\t%s\n",
            t,n,f,so,(f?100*so/f:0),(st?100*so/st:0),
            ao,(f?100*ao/f:0),(at?100*ao/at:0),
            fo,(f?100*fo/f:0),selected
        }' >> "$statistics_dir/coverage_by_gff_type.tsv"
    done

# Union-level accounting avoids double-counting bases shared by feature types.
awk -F '\t' 'BEGIN{OFS="\t"} $4!="region" {print $1,$2,$3}' "$tmp_dir/features.by_type.bed" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - \
  > "$tmp_dir/non_region_annotations.union.bed"
awk -F '\t' 'BEGIN{OFS="\t"} $4=="mobile_genetic_element" || tolower($4)~/repeat/ {print $1,$2,$3}' \
  "$tmp_dir/features.by_type.bed" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - \
  > "$tmp_dir/script_selected.union.bed"

printf 'mask_set\tmask_bp\toverlap_non_region_gff_bp\tpct_in_non_region_gff\toverlap_script_selected_bp\tpct_in_script_selected\n' \
  > "$statistics_dir/mask_set_summary.tsv"
for label in source final newly_added; do
  case "$label" in
    source) bed="$interval_dir/source_softmask.bed" ;;
    final) bed="$interval_dir/final_softmask.bed" ;;
    newly_added) bed="$interval_dir/newly_added_softmask.bed" ;;
  esac
  total="$(bed_bp "$bed")"
  all_overlap="$(overlap_bp "$bed" "$tmp_dir/non_region_annotations.union.bed")"
  selected_overlap="$(overlap_bp "$bed" "$tmp_dir/script_selected.union.bed")"
  awk -v label="$label" -v total="$total" -v all="$all_overlap" -v selected="$selected_overlap" \
    'BEGIN{OFS="\t"; printf "%s\t%d\t%d\t%.4f\t%d\t%.4f\n", label,total,all,(total?100*all/total:0),selected,(total?100*selected/total:0)}' \
    >> "$statistics_dir/mask_set_summary.tsv"
done

# Exact combinations of selected feature types provide non-double-counted
# attribution where repeat annotations overlap one another.
selected_types=(CENTROMERIC_REPEAT LONG_TERMINAL_REPEAT TELOMERIC_REPEAT mobile_genetic_element repeat_region)
selected_beds=()
for feature_type in "${selected_types[@]}"; do
  type_bed="$tmp_dir/selected.$feature_type.bed"
  awk -F '\t' -v t="$feature_type" 'BEGIN{OFS="\t"} $4==t{print $1,$2,$3}' \
    "$tmp_dir/features.by_type.bed" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > "$type_bed"
  selected_beds+=("$type_bed")
done
bedtools multiinter -i "${selected_beds[@]}" -names "${selected_types[@]}" \
  > "$tmp_dir/selected_type_combinations.bed"

printf 'mask_set\tfeature_type_combination\toverlap_bp\tpct_of_mask_set\n' \
  > "$statistics_dir/selected_type_combinations.tsv"
for label in source newly_added; do
  case "$label" in
    source) bed="$interval_dir/source_softmask.bed" ;;
    newly_added) bed="$interval_dir/newly_added_softmask.bed" ;;
  esac
  total="$(bed_bp "$bed")"
  bedtools intersect -a "$bed" -b "$tmp_dir/selected_type_combinations.bed" -wo \
    | awk -v label="$label" -v total="$total" 'BEGIN{OFS="\t"} {bp[$8]+=$NF} END{
        for (combination in bp)
          printf "%s\t%s\t%d\t%.4f\n",label,combination,bp[combination],100*bp[combination]/total
      }' \
    | LC_ALL=C sort -t $'\t' -k1,1 -k3,3nr \
    >> "$statistics_dir/selected_type_combinations.tsv"
done
