for k in {1..6}; do
	closestBed -d -a \
		../sd_meth/biser_TAIR10_rmmasked5/time-point."${k}".evolution.bed \
		-b ../structure/genebody.bed -t all \
		| awk '{
			if ($NF == 0) {
        if ($9 <= $2) {
            start = 0
        } else {
            start = $9 - $2
        }
        if ($10 >= $3) {
            end = $3 - $2
        } else {
            end = $10 - $2
        }
        printf "Chr%s(+):%d-%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, start, end, $4, $5, $6, $7, $11
			}
	}' >biser_TAIR10_rmmasked5/time-point."${k}".GB.list
	perl ../region_relative_pick.pl \
		-f biser_TAIR10_rmmasked5/time-point."${k}".muscle.fa -r biser_TAIR10_rmmasked5/time-point."${k}".GB.list \
		| awk '$3-$2>200' \
		| awk -F "," '$NF-$(NF-1)>200' \
		| sort -k1,1 -k2,2n >biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.bed
	closestBed -d -a \
		biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.bed \
		-b ../structure/genebody.bed -t all \
		| awk '$NF==0' \
		| perl ../bed_coverage.pl >biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.paired.bed
	perl -F'\t' -lane '($a, $b) = ((split(/\|/, $F[3]))[0], (split(/\|/, $F[7]))[0]); ($a, $b) = sort ($a, $b); print "$a\t$b"' \
		biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.paired.bed \
		| sort | uniq >biser_TAIR10_rmmasked5/time-point."${k}".paired.list

	while IFS=$'\t' read -r a b; do
		awk -va="$a.1" 'BEGIN {
				OFS = "\t"
				quote = "\047"
			}
			$1 == a {
				domains = (domains == "" ? quote $2 quote : domains ", " quote $2 quote)
			}
			END {
				if (domains != "") {
					print a, "[" domains "]"
				}
			}' ../data/gene_domains.tsv >tmp.1.tsv
		awk -va="$b.1" 'BEGIN {
				OFS = "\t"
				quote = "\047"
			}
			$1 == a {
				domains = (domains == "" ? quote $2 quote : domains ", " quote $2 quote)
			}
			END {
				if (domains != "") {
					print a, "[" domains "]"
				}
			}' ../data/gene_domains.tsv >tmp.2.tsv
		if [ -s tmp.1.tsv ] && [ -s tmp.2.tsv ]; then
			python3 ../domains_compare.py \
				-i tmp.1.tsv -I tmp.2.tsv \
				-s ../data/domain_score_matrix.tsv \
				-o biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.tmp
		fi
		rm tmp.1.tsv tmp.2.tsv
	done <biser_TAIR10_rmmasked5/time-point."${k}".paired.list

	awk -F'\t' '{sub(/\.1$/, "", $1); sub(/\.1$/, "", $2); print}' OFS='\t' \
		biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.tmp \
		>biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.tsv
	perl ../tsv_tiny1.pl \
		<biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.paired.bed \
		| sort | uniq \
		>biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.paired.list
	tsv-join biser_TAIR10_rmmasked5/time-point."${k}".muscle.GB.paired.list \
		-f biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.tsv \
		-k 1,2 -a 3 \
		>biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.score.tsv
done

for k in {1..6}; do
	awk '$NF>0' \
		biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.score.tsv \
		>biser_TAIR10_rmmasked5/time-point."${k}".paired.score.tsv
done
