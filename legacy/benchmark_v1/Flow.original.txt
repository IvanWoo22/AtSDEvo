### Prepare the reference sequences of TAIR10 and ColCEN and other info data.

```shell
mkdir "preref" && cd "preref" || exit
cp ../data/TAIR10_chr_all.fas.gz .
cp ../data/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz .
cp ../data/Col-CEN_v1.2.fasta.gz .

#Or, re-download:  
#wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
#wget https://github.com/schatzlab/Col-CEN/blob/main/v1.2/Col-CEN_v1.2.fasta.gz
#ARABIDOPSIS.ORG has been upgraded recently and might need to be downloaded manually.

pigz -dc TAIR10_chr_all.fas.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >TAIR10_unmasked.fa
rm Chr*.fa
pigz -dc Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz | sed 's/>/>Chr/g' | faops split-name stdin .
rm ChrPt.fa ChrMt.fa
cat Chr*.fa >TAIR10_E58masked.fa
rm Chr*.fa
pigz -dc Col-CEN_v1.2.fasta.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >ColCEN_unmasked.fa
rm Chr*.fa
rm TAIR10_chr_all.fas.gz Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz Col-CEN_v1.2.fasta.gz
```

### Mask the reference sequences.

```shell
for i in TAIR10; do
	mkdir ${i}_rm && cd ${i}_rm || exit
	faops split-name ../${i}_unmasked.fa .
	parallel -j 5 '
		RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
		RM2Bed.py -d . Chr{}.fa.out
		trf Chr{}.fa 2 5 7 80 10 40 500 -l 10 -h -ngs >Chr{}.fa.dat
		dustmasker -in Chr{}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr{}.fa_dust.bed
	' ::: {1..5}
	cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
	cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
	python ../../trf_merge.py Chr{1..5}.fa.dat trf.out.bed

	cat trf.out.bed dust.out.bed repeatmasker.out.bed \
		| cut -f 1-3 >>tmp.msk.bed
	cut -f 1-3 tmp.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -i - \
		| awk '$3-$2 > 2000 {print $0}' \
		| bedtools merge -d 100 -i - >tmp1.msk.bed
	cut -f 1-3 tmp.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -d 10 -i - >tmp2.msk.bed
	cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -i - \
		| seqtk seq -l 50 -M /dev/stdin ../${i}_unmasked.fa >../${i}_rmmasked.fa
	cd ..
done

perl ../../gff2bed_TE.pl ../../data/Araport11.cleanformat.gff >TE.bed
rm -r *_rm
```

### Masked coverage.

```shell
echo ",chr length,masked size,masked coverage" >masked_cover.csv
mkdir "smcover" && cd "smcover" || exit

egaz prepseq ../TAIR10_unmasked.fa -o .
for i in TAIR10_unmasked TAIR10_rmmasked TAIR10_rmmasked2 TAIR10_E58masked; do
	{
		echo -ne "${i},"
		faops masked ../${i}.fa | spanr cover stdin \
			| spanr stat --all chr.sizes stdin \
			| tail -1
	} >>../masked_cover.csv
done
rm chr.* Chr*
egaz prepseq ../ColCEN_unmasked.fa -o .
{
	awk '{sum+=$2};END{print "ColCEN_unmasked," sum ",0,0.0000"}' chr.sizes
	echo -ne "ColCEN_rmmasked,"
	faops masked ../ColCEN_rmmasked.fa | spanr cover stdin \
		| spanr stat --all chr.sizes stdin \
		| tail -1
} >>../masked_cover.csv
cd .. && rm -rf "smcover"
cd ..
```

`cat preref/masked_cover.csv`

|                      | chr length | masked size | masked coverage |
|:---------------------|-----------:|------------:|----------------:|
| TAIR10_unmasked      |  119146348 |      186207 |          0.0016 |
| TAIR10_rmmasked      |  119146348 |    21324693 |          0.1790 |
| TAIR10_rmmasked+dust |  119146348 |    26025948 |          0.2184 |
| TAIR10_E58masked     |  119146348 |    38123558 |          0.3200 |
| ColCEN_unmasked      |  131559676 |           0 |          0.0000 |
| ColCEN_rmmasked      |  131559676 |    33290733 |          0.2530 |
| ColCEN_rmmasked+dust |  131559676 |    38071449 |          0.2894 |

### `biser` for SD.

```shell
mkdir "biser" && cd "biser" || exit
parallel -j 6 '
	mkdir {} && cd {} || exit
	biser -t 4 -o biser_out ../../preref/{}.fa
	awk '\''{print $1"("$9"):"$2"-"$3"\t"$4"("$10"):"$5"-"$6}'\'' biser_out \
		| linkr sort stdin \
		| linkr clean stdin -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.95 \
		-o links.merge.tsv
	linkr clean links.sort.clean.tsv \
		-r links.merge.tsv --bundle 500 \
		-o links.clean.tsv
	linkr connect links.clean.tsv -r 0.05 \
		| linkr filter stdin -r 0.05 -o links.temp.tsv
	perl ../../check_overlap.pl --input links.temp.tsv --exclude --overlap 0.1 >links.tsv
	egaz prepseq ../../preref/{}.fa -o .
	perl -nla -F"\t" -e "print for @F" <links.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv \
			| perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " \
				>>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	rm -rf Chr*.fa chr* ./*.temp.yml ./links.merge.tsv ./links.sort.clean.tsv
' ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked
cd ..
```

### `lastz` for SD.

```shell
mkdir "lastz" && cd "lastz" || exit
echo 'strain,strain_id,species,species_id,genus,genus_id,family,family_id,order,order_id
Atha,3702,"Arabidopsis thaliana",3702,Arabidopsis,3701,Brassicaceae,3700,Brassicales,3699' >ensembl_taxon.csv
parallel -j 6 '
	mkdir {}
	egaz prepseq ../preref/{}.fa -o {}/
	egaz template {}/ --self -o {}/ --taxon ./ensembl_taxon.csv --circos --parallel 16 -v
	cd {} || exit
	mkdir -p Pairwise
	egaz lastz \
		--isself --set set01 -C 0 \
		--parallel 4 --verbose \
		. . \
		-o Pairwise
	egaz lpcnam \
		--parallel 4 --verbose \
		. . \
		Pairwise
	mkdir -p Results
	mkdir -p Processing
	ln -s "$(pwd)"/chr.fasta Processing/genome.fa
	cp -f "$(pwd)"/chr.sizes Processing/chr.sizes
	cd Processing || exit
	fasops axt2fas \
		../Pairwise/axtNet/*.axt.gz \
		-l 1000 -s chr.sizes -o stdout >axt.fas
	fasops separate axt.fas -o . --nodash -s .sep.fasta
	egaz exactmatch target.sep.fasta genome.fa --length 500 --discard 50 -o replace.target.tsv
	fasops replace axt.fas replace.target.tsv -o axt.target.fas
	egaz exactmatch query.sep.fasta genome.fa --length 500 --discard 50 -o replace.query.tsv
	fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas
	fasops covers axt.correct.fas -o axt.correct.yml
	spanr split axt.correct.yml -s .temp.yml -o .
	spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
	spanr stat chr.sizes axt.union.yml -o union.csv
	fasops links axt.correct.fas -o stdout | perl -nl -e "s/(target|query)\.//g; print;" >links.lastz.tsv
	fasops separate axt.correct.fas --nodash --rc -o stdout \
		| perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" \
		| faops filter -u stdin stdout \
		| faops filter -n 250 stdin stdout \
			>axt.gl.fasta
	egaz blastn axt.gl.fasta genome.fa -o axt.bg.blast --parallel 8
	egaz blastmatch axt.bg.blast -c 0.95 -o axt.bg.region --parallel 8
	samtools faidx genome.fa -r axt.bg.region --continue \
		| perl -p -e "/^>/ and s/:/(+):/" >axt.bg.fasta
	cat axt.gl.fasta axt.bg.fasta | faops filter -u stdin stdout \
		| faops filter -n 250 stdin stdout >axt.all.fasta
	egaz blastn axt.all.fasta axt.all.fasta -o axt.all.blast --parallel 8
	egaz blastlink axt.all.blast -c 0.95 -o links.blast.tsv --parallel 8
	linkr sort links.lastz.tsv links.blast.tsv -o links.sort.tsv
	linkr clean links.sort.tsv -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.9 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 500 -o links.clean.tsv
	linkr connect links.clean.tsv -r 0.8 -o links.connect.tsv
	linkr filter links.connect.tsv -r 0.8 -o links.temp.tsv
	perl ../../../check_overlap.pl --input links.temp.tsv --exclude --overlap 0.15 >links.filter.tsv
	fasops create links.filter.tsv -g genome.fa -o multi.temp.fas
	fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o multi.refine.fas
	fasops links multi.refine.fas -o stdout | linkr sort stdin -o stdout | linkr filter stdin -n 2-50 -o links.refine.tsv
	fasops links multi.refine.fas -o stdout --best | linkr sort stdin -o links.best.tsv
	fasops create links.best.tsv -g genome.fa --name Atha -o pair.temp.fas
	fasops refine pair.temp.fas --msa mafft -p 16 -o pair.refine.fas
	perl -nla -F"\t" -e "print for @F" <links.refine.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.refine.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv \
			| perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " \
				>>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	cp cover.yml ../cover.yml
	cp copy.yml ../copy.yml
	mv cover.yml.csv ../cover.csv
	mv copy.csv ../copy.csv
	cp links.refine.tsv ../links.tsv
	mv multi.refine.fas ../multi.fas
	mv pair.refine.fas ../pair.fas
	cd ..
	rm -rf Chr*.fa chr* *.sh Processing Pairwise Results
' ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked
cd ..
```

### Detect functional structure of genes.

```shell
mkdir structure
perl gff2bed_split.pl \
	data/Araport11_gene_type.txt \
	data/Araport11_cleanformat.gff \
	structure

for i in promoter genebody gene; do
	sort -k1,1 -k2,2n structure/${i}.raw.bed >structure/${i}.bed
	for j in pseudogene noncoding_rna novel_transcribed_region protein_coding; do
		grep ${j} structure/${i}.bed >structure/${j}.${i}.bed
	done
done

#awk '$1~/^[0-9]+$/{print $1 "\t" $2 "\t" $2 "\t" $7 "\t" $5 "\t" $3}' \
#	~/fat/NJU_seq/mrna.scored.Nm.{Pst,MgCl2}.tsv \
#	| sort -k1,1 -k2,2n \
#	| uniq >data/Atha.mrna.Nm.bed
cut -f 4 data/Atha.mrna.Nm.bed \
	| awk -F "/" '{print $1}' \
	| sort | uniq >data/Atha.mrna.Nm.genelist.tsv
grep -f Atha.mrna.Nm.gene.list structure/gene.bed >structure/withNm.gene.bed
rm structure/*.raw.bed
```

```shell
mkdir "sd_meth" && cd "sd_meth" || exit
for i in Alyrata Ahalleri Bstricta Rislandica BrapaFPsc Esyriacum Dstrictus Cpapaya Cviolacea Tcacao Graimondii; do
	ln -s "$PWD"/../phyt/Atha_${i}.yml Atha_${i}.yml
done

parallel -j 6 '
	mkdir {1}_{2}
	linkr filter ../{1}/{2}/Atha.links.tsv -n 2 -o {1}_{2}/links.2copy.tmp
	perl ../split_lines.pl {1}_{2}/links.2copy.tmp {1}_{2}/links.2copy.tsv
	rm {1}_{2}/links.2copy.tmp
	perl ../classify_links.pl {1}_{2}/links.2copy.tsv \
		Atha_{Alyrata,Ahalleri,Bstricta,Rislandica,BrapaFPsc,Esyriacum,Dstrictus,Cpapaya,Cviolacea,Tcacao,Graimondii}.yml \
		0.5 >{1}_{2}/time-point.raw.tsv
	cut -f 1,2 {1}_{2}/time-point.raw.tsv >{1}_{2}/time-point.0.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 3,4 >{1}_{2}/time-point.1.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 5,6 >{1}_{2}/time-point.2.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 7,8,9 >{1}_{2}/time-point.3.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 10,11 >{1}_{2}/time-point.4.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 12,13 >{1}_{2}/time-point.5.tmp
	paste {1}_{2}/time-point.{0..5}.tmp >{1}_{2}/time-point.tsv
	rm {1}_{2}/time-point.{0..5}.tmp {1}_{2}/time-point.raw.tsv
' ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked
```

```shell
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 0 "\t" $1}' >time1-0.tmp
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time1-1.tmp
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time1-2.tmp
		cat time1-{0,1,2}.tmp >${i}_${j}/time1.txt
		for k in {2..5}; do
			cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
				| awk -va="$k" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (k = (a+1); k <= NF; k++) {
						if ($k == 0) {
								count0++;
						} else if ($k == 1) {
								count1++;
						} else if ($k == 2) {
								count2++;
						} else if ($k == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 == 0 && count2 == 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$k" '$a==3{print $2 $3 $4 $5 $6 "\t" 0 "\t" $1}' >time"${k}"-0.tmp
			cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
				| awk -va="$k" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (k = (a+1); k <= NF; k++) {
						if ($k == 0) {
								count0++;
						} else if ($k == 1) {
								count1++;
						} else if ($k == 2) {
								count2++;
						} else if ($k == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$k" '$a==3{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time"${k}"-1.tmp
			cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
				| awk -va="$k" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (k = (a+1); k <= NF; k++) {
						if ($k == 0) {
								count0++;
						} else if ($k == 1) {
								count1++;
						} else if ($k == 2) {
								count2++;
						} else if ($k == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$k" '$a==3{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time"${k}"-2.tmp
			cat time"${k}"-{0,1,2}.tmp >${i}_${j}/time"${k}".txt
		done
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk '$NF==3{print $2 $3 $4 $5 $6 $7 "\t" 0 "\t" $1}' >${i}_${j}/time6.txt
		rm time*.tmp
		awk -va=${i} -vb=${j} '{print b "\t" a "\t" ($3-$2) "\n" b "\t" a "\t" ($6-$5)}' \
			${i}_${j}/links.2copy.tsv >>links2_length_distribution.tsv
		for k in {1..6}; do
			awk -va=${i} -vb=${j} -vc="${k}" '{sum+=$3};END{print b "\t" a "\t" c "\t" sum}' ${i}_${j}/time"${k}".txt >>links2_timedistribution.tsv
		done
		awk -F"\t" -va=${i} -vb=${j} '$1==b&&$2==a{sum+=$4};END{print b "\t" a "\t" sum}' links2_timedistribution.tsv >>links2_timedistributionallcount.tsv
		wc -l <${i}_${j}/links.2copy.tsv >>links2_allcount.tsv
	done
done
paste links2_timedistributionallcount.tsv links2_allcount.tsv >links2_counts_in_timeline.tsv
rm links2_timedistributionallcount.tsv links2_allcount.tsv
```

```shell
parallel -j 24 "
	perl ../fetch_time_point.pl {1}_{2}/time{3}.txt {1}_{2}/time-point.tsv \\
		>{1}_{2}/time-point.{3}.tsv
	awk '\$7==1{print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" NR \"\\t1\\n\" \$4 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" NR \"\\t2\"};
		\$7==0{print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" NR \"\\t0\\n\" \$4 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" NR \"\\t0\"}' \\
		{1}_{2}/time-point.{3}.tsv \\
		| sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.sort.bed
	awk '{print \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" NR \"\\t\" \$8 \"\\n\" \$4 \"\\t\" \$5 \"\\t\" \$6 \"\\t\" NR \"\\t\" \$8}' \\
		{1}_{2}/time-point.{3}.tsv \\
		| sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.evolution.bed
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}
```

```shell
parallel -j 24 "
	perl ../arg_meth_link.pl \\
		../../MASED/Memory/AT.beta.1.tsv \\
		{1}_{2}/time-point.{3}.sort.bed \\
		{1}_{2}/time-point.{3}.beta.bed
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
				${i}_${j}/time-point.${k}.beta.bed >>${i}_${j}/time-point.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
				${i}_${j}/time-point.${k}.beta.bed >>time-point.beta.all.bed
		done
	done
done
```

```shell
parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.prom.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.prom.bed \
		{1}_{2}/time-point.{3}.prom.beta.bed 6
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.prom.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.prom.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point."${k}".prom.beta.bed \
				>>${i}_${j}/time-point.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point."${k}".prom.beta.bed \
				>>time-point.prom.beta.all.bed
		done
	done
done
```

```shell
parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/pseudogene.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.pseudo.prom.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.pseudo.prom.bed \
		{1}_{2}/time-point.{3}.pseudo.prom.beta.bed 6
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.pseudo.prom.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.pseudo.prom.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point."${k}".pseudo.prom.beta.bed \
				>>${i}_${j}/time-point.pseudo.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point."${k}".pseudo.prom.beta.bed \
				>>time-point.pseudo.prom.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/protein_coding.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.pc.prom.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.pc.prom.bed \
		{1}_{2}/time-point.{3}.pc.prom.beta.bed 6
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.pc.prom.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.pc.prom.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pc.prom.beta.bed \
				>>${i}_${j}/time-point.pc.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pc.prom.beta.bed \
				>>time-point.pc.prom.beta.all.bed
		done
	done
done
```

```shell
rm links2_Nm_cover_in_timeline.tsv
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		for k in {1..6}; do
			awk '{print $1 ":" $2 "-" $3}' ${i}_${j}/time-point."${k}".sort.bed \
				| spanr cover stdin -o temp.yml
			echo -ne "${j}\t${i}\t${k}\t" >>links2_Nm_cover_in_timeline.tsv
			spanr stat ../../MASED/revis/RMasked/Ensembl/Atha/chr.sizes \
				--all temp.yml | awk -F "," 'NR==2{printf $2 "\t"}' >>links2_Nm_cover_in_timeline.tsv
			rm temp.yml
			closestBed -d -a ${i}_${j}/time-point."${k}".sort.bed \
				-b ../Atha.mrna.Nm.bed -t all \
				| awk '$NF==0' | cut -f 6-8 \
				| sort | uniq | wc -l >>links2_Nm_cover_in_timeline.tsv
		done
	done
done
```