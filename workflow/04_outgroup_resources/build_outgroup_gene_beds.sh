#!/usr/bin/env bash
set -euo pipefail

inventory_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="$inventory_dir/selected_outgroups.tsv"
out_dir="$inventory_dir/gene_coordinates_bed"
mkdir -p "$out_dir"

awk -F '\t' 'NR>1{print $2"\t"$8}' "$manifest" |
while IFS=$'\t' read -r code coordinate_gff; do
  # Historical MCScanX coordinate files are four-column, 1-based records:
  # sequence, transcript ID, endpoint 1, endpoint 2. Normalize orientation and
  # emit standard 0-based half-open BED6 without changing sequence names.
  awk -F '\t' 'BEGIN{OFS="\t"} NF>=4 {
      if($3 <= $4){start=$3-1; end=$4; strand="+"}
      else{start=$4-1; end=$3; strand="-"}
      print $1,start,end,$2,0,strand
    }' "$coordinate_gff" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    > "$out_dir/$code.genes.bed"
done
