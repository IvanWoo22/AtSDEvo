#!/usr/bin/env bash
set -euo pipefail

reference_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src_dir="$reference_dir/source_archive/dataset/ncbi_dataset/data/GCA_978657495.1"
archive_genome="$src_dir/GCA_978657495.1_TAIR12_genomic.fna"
archive_gff="$src_dir/genomic.gff"
out_dir="$reference_dir/prepared_data"
meta_dir="$reference_dir/qc_reports"
genome="$out_dir/TAIR12.Col-CC.source_softmasked.fa"
gff="$out_dir/TAIR12.Col-CC.annotation.gff3"

mkdir -p "$out_dir" "$meta_dir"

# Normalize the five INSDC accessions to the chromosome names explicitly
# recorded in sequence_report.jsonl. The archive files remain untouched.
awk '
  BEGIN {
    map["OZ408683.1"]="Chr1"; map["OZ408684.1"]="Chr2";
    map["OZ408685.1"]="Chr3"; map["OZ408686.1"]="Chr4";
    map["OZ408687.1"]="Chr5"
  }
  /^>/ {
    name=substr($1,2); if(name in map) $1=">" map[name]; print; next
  }
  {print}
' "$archive_genome" > "$genome"

awk '
  BEGIN {
    map["OZ408683.1"]="Chr1"; map["OZ408684.1"]="Chr2";
    map["OZ408685.1"]="Chr3"; map["OZ408686.1"]="Chr4";
    map["OZ408687.1"]="Chr5"
  }
  {
    for(old in map) gsub(old,map[old]); print
  }
' "$archive_gff" > "$gff"

samtools faidx "$genome"

# GFF3 is 1-based closed; BED is 0-based half-open. Keep explicit repeat and
# mobile-element features. Existing lowercase sequence is retained by
# bedtools maskfasta, so this supplements rather than erases the archive mask.
awk -F '\t' 'BEGIN{OFS="\t"}
  !/^#/ && ($3 == "mobile_genetic_element" || tolower($3) ~ /repeat/) {
    print $1, $4 - 1, $5
  }' "$gff" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > "$out_dir/TAIR12.annotated_repeats.merged.bed"

bedtools maskfasta -soft \
  -fi "$genome" \
  -bed "$out_dir/TAIR12.annotated_repeats.merged.bed" \
  -fo "$out_dir/TAIR12.Col-CC.annotation_softmasked.fa"

samtools faidx "$out_dir/TAIR12.Col-CC.annotation_softmasked.fa"

md5sum "$archive_genome" "$archive_gff" "$genome" "$gff" \
  "$out_dir/TAIR12.annotated_repeats.merged.bed" \
  "$out_dir/TAIR12.Col-CC.annotation_softmasked.fa" \
  > "$meta_dir/md5sums.txt"

perl -ne '
  next if /^>/;
  $upper += tr/ACGT//;
  $lower += tr/acgt//;
  $n += tr/Nn//;
  END { print "uppercase_ACGT\t$upper\nlowercase_acgt\t$lower\nN_or_n\t$n\n" }
' "$genome" > "$meta_dir/source_fasta_case_counts.tsv"

perl -ne '
  next if /^>/;
  $upper += tr/ACGT//;
  $lower += tr/acgt//;
  $n += tr/Nn//;
  END { print "uppercase_ACGT\t$upper\nlowercase_acgt\t$lower\nN_or_n\t$n\n" }
' "$out_dir/TAIR12.Col-CC.annotation_softmasked.fa" \
  > "$meta_dir/prepared_fasta_case_counts.tsv"

awk 'BEGIN{n=0; bp=0} {n++; bp += $3-$2} END{
  print "merged_intervals\t" n;
  print "annotated_repeat_bp\t" bp
}' "$out_dir/TAIR12.annotated_repeats.merged.bed" \
  > "$meta_dir/annotated_repeat_bed_stats.tsv"
