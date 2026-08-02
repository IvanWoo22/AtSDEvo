mkdir "preref" && cd "preref" || exit
cp ../data/TAIR10_chr_all.fas.gz .
cp ../data/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz .
pigz -dc TAIR10_chr_all.fas.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >TAIR10_unmasked.fa
rm Chr*.fa
pigz -dc Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz | sed 's/>/>Chr/g' | faops split-name stdin .
rm ChrPt.fa ChrMt.fa
cat Chr*.fa >TAIR10_E58masked.fa
rm Chr*.fa
perl ../gff2bed_TE.pl ../data/Araport11.cleanformat.gff >TE.bed

mkdir TAIR10_rm1 && cd TAIR10_rm1 || exit
faops split-name ../TAIR10_unmasked.fa .
parallel -j 5 '
	RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
	RM2Bed.py -d . Chr{}.fa.out
	' ::: {1..5}
cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
cut -f 1-3 repeatmasker.out.bed >>tmp.msk.bed
bedtools sort -i tmp.msk.bed \
	| awk '$3-$2>2000{print $0}' \
	| bedtools merge -d 100 -i - >tmp1.msk.bed
bedtools sort -i tmp.msk.bed \
	| bedtools merge -d 10 -i - >tmp2.msk.bed
cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked1.fa
cd ..

mkdir TAIR10_rm2 && cd TAIR10_rm2 || exit
faops split-name ../TAIR10_unmasked.fa .
parallel -j 5 '
	RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
	RM2Bed.py -d . Chr{}.fa.out
	dustmasker -in Chr{}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr{}.fa_dust.bed
	' ::: {1..5}
cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
cat dust.out.bed repeatmasker.out.bed \
	| cut -f 1-3 >>tmp.msk.bed
bedtools sort -i tmp.msk.bed \
	| awk '$3-$2>2000{print $0}' \
	| bedtools merge -d 100 -i - >tmp1.msk.bed
bedtools sort -i tmp.msk.bed \
	| bedtools merge -d 10 -i - >tmp2.msk.bed
cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked2.fa
cd ..

mkdir TAIR10_rm3 && cd TAIR10_rm3 || exit
faops split-name ../TAIR10_unmasked.fa .
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
bedtools sort -i tmp.msk.bed \
	| awk '$3-$2>2000{print $0}' \
	| bedtools merge -d 100 -i - >tmp1.msk.bed
bedtools sort -i tmp.msk.bed \
	| bedtools merge -d 10 -i - >tmp2.msk.bed
cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked3.fa
cd ..

mkdir TAIR10_rm4 && cd TAIR10_rm4 || exit
faops split-name ../TAIR10_unmasked.fa .
parallel -j 5 '
	RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
	RM2Bed.py -d . Chr{}.fa.out
	dustmasker -in Chr{}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr{}.fa_dust.bed
	' ::: {1..5}
cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
cat dust.out.bed repeatmasker.out.bed ../TE.bed \
	| cut -f 1-3 >tmp.msk.bed
bedtools sort -i tmp.msk.bed \
	| awk '$3-$2>2000{print $0}' \
	| bedtools merge -d 100 -i - >tmp1.msk.bed
bedtools sort -i tmp.msk.bed \
	| bedtools merge -d 10 -i - >tmp2.msk.bed
cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked4.fa
cd ..

mkdir TAIR10_rm5 && cd TAIR10_rm5 || exit
faops split-name ../TAIR10_unmasked.fa .
parallel -j 5 '
	RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
	RM2Bed.py -d . Chr{}.fa.out
	trf Chr{}.fa 2 5 7 80 10 40 500 -l 10 -h -ngs >Chr{}.fa.dat
	dustmasker -in Chr{}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr{}.fa_dust.bed
	' ::: {1..5}
cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
python ../../trf_merge.py Chr{1..5}.fa.dat trf.out.bed
cat trf.out.bed dust.out.bed repeatmasker.out.bed ../TE.bed \
	| cut -f 1-3 >tmp.msk.bed
bedtools sort -i tmp.msk.bed \
	| awk '$3-$2>2000{print $0}' \
	| bedtools merge -d 100 -i - >tmp1.msk.bed
bedtools sort -i tmp.msk.bed \
	| bedtools merge -d 10 -i - >tmp2.msk.bed
cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked5.fa
cd ..

echo ",chr length,masked size,masked coverage" >masked_cover.csv
mkdir "smcover" && cd "smcover" || exit
egaz prepseq ../TAIR10_unmasked.fa -o .
for i in TAIR10_rmmasked{1..5}; do
	{
		echo -ne "${i},"
		faops masked ../"${i}".fa | spanr cover stdin \
			| spanr stat --all chr.sizes stdin \
			| tail -1
	} >>../masked_cover.csv
done
cd .. && rm -rf "smcover"
cd ..

mkdir "biser" && cd "biser" || exit
parallel -j 3 '
	mkdir {} && cd {} || exit
	biser -t 6 -o biser_out ../../preref/{}.fa
	awk '\''{print $1"("$9"):"$2"-"$3"\t"$4"("$10"):"$5"-"$6}'\'' biser_out | linkr sort stdin | linkr clean stdin -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.9 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 0 -o links.clean.tsv
	linkr connect links.clean.tsv -r 0.05 | linkr filter stdin -r 0.05 -o links.temp.tsv
	perl ../../check_overlap.pl --input links.temp.tsv --exclude --overlap 0.2 >links.tsv
	singularity run $HOME/egaz_master.sif egaz prepseq ../../preref/{}.fa -o .
	perl -nla -F"\t" -e "print for @F" <links.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv | perl -nl -e "
				@fields = grep {/\S+/} split /\s+/;
				next unless @fields == 2;
				next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
				printf qq{%s,%s\n}, \$1, \$fields[0];
			" >>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	rm -rf Chr*.fa chr* ./*.temp.*
' ::: TAIR10_rmmasked{1..5}
cd ..

mkdir "lastz" && cd "lastz" || exit
parallel -j 3 '
	mkdir {}
	singularity run $HOME/egaz_master.sif egaz prepseq ../preref/{}.fa -o {}/
	cd {} || exit
	mkdir -p Pairwise
	singularity run $HOME/egaz_master.sif egaz lastz --isself --set set01 -C 0 --parallel 4 --verbose . . -o Pairwise
	singularity run $HOME/egaz_master.sif egaz lnam --parallel 4 --verbose . . Pairwise
	mkdir -p Results
	mkdir -p Processing
	ln -s "$(pwd)"/chr.fasta Processing/genome.fa
	cd Processing || exit
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops axt2fas ../Pairwise/axtNet/*.axt.gz -l 1000 -s ../chr.sizes -o stdout >axt.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops separate axt.fas -o . --nodash -s .sep.fasta
	singularity run $HOME/egaz_master.sif egaz exactmatch target.sep.fasta genome.fa --length 500 --discard 50 -o replace.target.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops replace axt.fas replace.target.tsv -o axt.target.fas
	singularity run $HOME/egaz_master.sif egaz exactmatch query.sep.fasta genome.fa --length 500 --discard 50 -o replace.query.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops covers axt.correct.fas -o axt.correct.yml
	spanr split axt.correct.yml -s .temp.yml -o .
	spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
	spanr stat ../chr.sizes axt.union.yml -o ../union.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops links axt.correct.fas -o stdout | perl -nl -e "s/(target|query)\.//g; print;" >links.lastz.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops separate axt.correct.fas --nodash --rc -o stdout | perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" | faops filter -u stdin stdout | faops filter -n 250 stdin stdout >axt.gl.fasta
	singularity run $HOME/egaz_master.sif egaz blastn axt.gl.fasta genome.fa -o axt.bg.blast --parallel 8
	singularity run $HOME/egaz_master.sif egaz blastmatch axt.bg.blast -c 0.95 -o axt.bg.region --parallel 8
	samtools faidx genome.fa -r axt.bg.region --continue | perl -p -e "/^>/ and s/:/(+):/" >axt.bg.fasta
	cat axt.gl.fasta axt.bg.fasta | faops filter -u stdin stdout | faops filter -n 250 stdin stdout >axt.all.fasta
	singularity run $HOME/egaz_master.sif egaz blastn axt.all.fasta axt.all.fasta -o axt.all.blast --parallel 8
	singularity run $HOME/egaz_master.sif egaz blastlink axt.all.blast -c 0.95 -o links.blast.tsv --parallel 8
	linkr sort links.lastz.tsv links.blast.tsv -o links.sort.tsv
	linkr clean links.sort.tsv -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.9 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 0 -o links.clean.tsv
	linkr connect links.clean.tsv -r 0.9 -o links.connect.tsv
	linkr filter links.connect.tsv -r 0.9 -o links.filter.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops create links.filter.tsv -g genome.fa -o multi.temp.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o multi.refine.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops links multi.refine.fas -o stdout | linkr sort stdin -o stdout | linkr filter stdin -n 2-50 -o links.temp.tsv
	perl ../../../check_overlap.pl --input links.temp.tsv --exclude --overlap 0.05 >../links.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops create ../links.tsv -g genome.fa -o multi.temp.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o ../multi.fas
	cd ..
	perl -nla -F"\t" -e "print for @F" <links.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv | perl -nl -e "
				@fields = grep {/\S+/} split /\s+/;
				next unless @fields == 2;
				next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
				printf qq{%s,%s\n}, \$1, \$fields[0];
			" >>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.csv
	rm -rf Chr*.fa chr* *.temp.* Processing Pairwise Results
' ::: TAIR10_rmmasked{1..5}
cd ..

mkdir "sd_meth" && cd "sd_meth" || exit
mkdir "colink"
for i in Alyrata Ahalleri Bstricta Rislandica BrapaFPsc Esyriacum Dstrictus Cpapaya Cviolacea Tcacao Graimondii; do
	ln -s "$PWD"/../phyt/Atha_${i}.yml colink/Atha_${i}.yml
done
parallel -j 6 '
	mkdir {1}_{2}
	linkr filter ../{1}/{2}/links.tsv -n 2 -o {1}_{2}/links.2copy.tmp
	perl ../split_lines.pl {1}_{2}/links.2copy.tmp {1}_{2}/links.2copy.tsv
	rm {1}_{2}/links.2copy.tmp
	perl ../classify_links.pl {1}_{2}/links.2copy.tsv \
		colink/Atha_{Alyrata,Ahalleri,Bstricta,Rislandica,BrapaFPsc,Esyriacum,Dstrictus,Cpapaya,Cviolacea,Tcacao,Graimondii}.yml 0.5 \
		>{1}_{2}/time-point.raw.tsv
	cut -f 1,2,14 {1}_{2}/time-point.raw.tsv >{1}_{2}/time-point.0.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 3,4 >{1}_{2}/time-point.1.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 5,6 >{1}_{2}/time-point.2.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 7,8,9 >{1}_{2}/time-point.3.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 10,11 >{1}_{2}/time-point.4.tmp
	python ../tsv_summarize.py {1}_{2}/time-point.raw.tsv 12,13 >{1}_{2}/time-point.5.tmp
	paste {1}_{2}/time-point.{0..5}.tmp >{1}_{2}/time-point.tsv
	rm {1}_{2}/time-point.{0..5}.tmp {1}_{2}/time-point.raw.tsv
' ::: lastz biser ::: TAIR10_rmmasked{1..5}

for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		bash ../colink2time_point.sh "${i}" "${j}"
	done
done
paste links2_timedistributionallcount.tsv links2_allcount.tsv >links2_counts_in_timeline.tsv
rm links2_timedistributionallcount.tsv links2_allcount.tsv

parallel -j 24 '
	perl ../fetch_time_point.pl {1}_{2}/time{3}.txt {1}_{2}/time-point.tsv >{1}_{2}/time-point.{3}.tsv
	awk '\''$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2"}
					$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0"}'\'' \
		{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.sort.bed
	awk '\''$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2\t" $7 "\t" $9}
					$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0\t" $7 "\t" $9}'\'' \
		{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.evolution.bed
	perl ../arg_meth_link.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.sort.bed \
		{1}_{2}/time-point.{3}.beta.bed
' ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}

rm time-point.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
				${i}_"${j}"/time-point."${k}".beta.bed >>${i}_"${j}"/time-point.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
				${i}_"${j}"/time-point."${k}".beta.bed >>time-point.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.promoter.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.promoter.bed \
		{1}_{2}/time-point.{3}.promoter.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}

rm time-point.promoter.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.promoter.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".promoter.beta.bed \
				>>${i}_"${j}"/time-point.promoter.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".promoter.beta.bed \
				>>time-point.promoter.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/pseudogene.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.pseudo.promoter.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.pseudo.promoter.bed \
		{1}_{2}/time-point.{3}.pseudo.promoter.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}
rm time-point.pseudo.promoter.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.pseudo.promoter.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".pseudo.promoter.beta.bed \
				>>${i}_"${j}"/time-point.pseudo.promoter.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".pseudo.promoter.beta.bed \
				>>time-point.pseudo.promoter.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../promotor_intsec.pl \
		../structure/protein_coding.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}..promoter.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}..promoter.bed \
		{1}_{2}/time-point.{3}..promoter.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}
rm time-point..promoter.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point..promoter.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}"..promoter.beta.bed \
				>>${i}_"${j}"/time-point..promoter.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}"..promoter.beta.bed \
				>>time-point..promoter.beta.all.bed
		done
	done
done
rm links2_cover_in_timeline.tsv 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		for k in {1..6}; do
			awk '{print $1 ":" $2 "-" $3}' ${i}_"${j}"/time-point."${k}".sort.bed \
				| spanr cover stdin -o temp.yml
			echo -ne "${j}\t${i}\t${k}\t" >>links2_cover_in_timeline.tsv
			spanr stat ../../MASED/revis/RMasked/Ensembl/Atha/chr.sizes \
				--all temp.yml | awk -F "," 'NR==2{printf $2 "\t"}' >>links2_cover_in_timeline.tsv
			rm temp.yml
			closestBed -d -a ${i}_"${j}"/time-point."${k}".sort.bed \
				-b ../data/Atha.mrna.NmNeo.bed -t all \
				| awk '$NF==0' | cut -f 6-8 \
				| sort | uniq | wc -l >>links2_cover_in_timeline.tsv
		done
	done
done

sed 's/Chr//g' ../preref/TE.bed | sort -k1,1 -k2,2n >./TE.bed
parallel -j 20 "
	perl ../bed_intsec.pl \
		--diff --ref_bed ./TE.bed \
		--in_bed {1}_{2}/time-point.{3}.sort.bed \
		--output {1}_{2}/time-point.{3}.nonTE.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.nonTE.bed \
		{1}_{2}/time-point.{3}.nonTE.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}
rm time-point.nonTE.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.nonTE.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".nonTE.beta.bed \
				>>${i}_"${j}"/time-point.nonTE.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".nonTE.beta.bed \
				>>time-point.nonTE.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../bed_intsec.pl \
		--ref_bed ./TE.bed \
		--in_bed {1}_{2}/time-point.{3}.sort.bed \
		--output {1}_{2}/time-point.{3}.TE.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.TE.bed \
		{1}_{2}/time-point.{3}.TE.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}
rm time-point.TE.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.TE.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".TE.beta.bed \
				>>${i}_"${j}"/time-point.TE.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".TE.beta.bed \
				>>time-point.TE.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../bed_intsec.pl \
		--ref_bed ../structure/promoter.bed \
		--in_bed {1}_{2}/time-point.{3}.nonTE.bed \
		--output {1}_{2}/time-point.{3}.nonTE.promoter.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.nonTE.promoter.bed \
		{1}_{2}/time-point.{3}.nonTE.promoter.beta.bed 6
" ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}
rm time-point.nonTE.promoter.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		rm ${i}_"${j}"/time-point.nonTE.promoter.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".nonTE.promoter.beta.bed \
				>>${i}_"${j}"/time-point.nonTE.promoter.beta.bed
			awk -va=${i} -vb="${j}" -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_"${j}"/time-point."${k}".nonTE.promoter.beta.bed \
				>>time-point.nonTE.promoter.beta.all.bed
		done
	done
done

for i in {1..6}; do
	{
		wc -l <biser_TAIR10_rmmasked5/time-point."${i}".sort.bed
		closestBed -d -a biser_TAIR10_rmmasked5/time-point."${i}".sort.bed -b TE.bed -t all | awk '$NF==0' | cut -f 1-5 | sort | uniq | wc -l
		awk '$5==1' biser_TAIR10_rmmasked5/time-point."${i}".sort.bed | wc -l
		closestBed -d -a biser_TAIR10_rmmasked5/time-point."${i}".sort.bed -b TE.bed -t all | awk '$NF==0' | cut -f 1-5 | sort | uniq | awk '$5==1' | wc -l
		awk '$5==2' biser_TAIR10_rmmasked5/time-point."${i}".sort.bed | wc -l
		closestBed -d -a biser_TAIR10_rmmasked5/time-point."${i}".sort.bed -b TE.bed -t all | awk '$NF==0' | cut -f 1-5 | sort | uniq | awk '$5==2' | wc -l
		awk '$5==0' biser_TAIR10_rmmasked5/time-point."${i}".sort.bed | wc -l
		closestBed -d -a biser_TAIR10_rmmasked5/time-point."${i}".sort.bed -b TE.bed -t all | awk '$NF==0' | cut -f 1-5 | sort | uniq | awk '$5==0' | wc -l
	} >temp"${i}".tsv
done

cd ..
mkdir "sd_simi" && cd "sd_simi" || exit
for i in lastz biser; do
	for j in TAIR10_rmmasked{1..5}; do
		mkdir "${i}"_"${j}"
		for k in {1..6}; do
			mkdir tmp
			perl ../bed2links.pl \
				../sd_meth/"${i}"_"${j}"/time-point."${k}".evolution.bed \
				../preref/"${j}".fa \
				"${i}"_"${j}"/time-point."${k}".info.tsv tmp
			for fa_file in tmp/*.fa; do
				fasops refine "${fa_file}" --msa muscle -o "${fa_file}".out
			done
			cat tmp/*.fa.out >"${i}"_"${j}"/time-point."${k}".muscle.fa
			rm -r tmp
		done
	done
done

for k in {1..6}; do
	closestBed -d -a \
		../sd_meth/biser_TAIR10_rmmasked5/time-point."${k}".evolution.bed \
		-b ../sd_meth/TE.bed -t all \
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
	}' >biser_TAIR10_rmmasked5/time-point."${k}".TE.list
	perl ../region_blank_pick.pl \
		-f biser_TAIR10_rmmasked5/time-point."${k}".muscle.fa \
		-r biser_TAIR10_rmmasked5/time-point."${k}".TE.list \
		>biser_TAIR10_rmmasked5/time-point."${k}".muscle.TE.tsv
	rm biser_TAIR10_rmmasked5/time-point."${k}".TE.list
done

for k in {1..6}; do
	closestBed -d -a \
		../sd_meth/biser_TAIR10_rmmasked5/time-point."${k}".evolution.bed \
		-b ../structure/protein_coding.genebody.bed -t all \
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
		-b ../structure/protein_coding.genebody.bed -t all \
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
	rm .GB.list biser_TAIR10_rmmasked5/time-point."${k}".paired.compared.tmp
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

awk '{print $0 "\t" $2 "-" $3}' ../structure/promoter.bed >promoter.bed

perl ../arg_meth_link_neo.pl \
	../../MASED/Memory/AT.beta.1.tsv \
	promoter.bed \
	promoter.beta.bed 5
