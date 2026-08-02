awk '$3 == "CDS" {print $9 ";"}' genomic.gff | head -2
awk '$3 == "mRNA" {print $9 ";"}' genomic.gff | head -2
grep "^>" protein.faa | head -10000 | tail -2
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' genomic.gff >gene.gff
head -2 gene.gff

awk '$3 == "CDS" {print $9 ";"}' genomic.gff \
	| perl -ne '/ID=cds-([^;]+);*Parent=rna-gnl\|JCVI\|mRNA\.([^;]+)/;print"$2\t$1\n";' \
	| uniq >cds.list
awk '$3 == "mRNA" {print $9 ";"}' genomic.gff \
	| perl -ne '/ID=rna-gnl\|JCVI\|mRNA\.([^;]+);*Parent=gene-([^;]+)/;print"$1\t$2\n";' \
	| uniq >mrna.list
tsv-join cds.list -f mrna.list -k 1 -a 2 | cut -f 2,3 >alias.txt
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' genomic.gff >gene.gff
perl ../../gff_pep.pl \
	--abbr "Csat" --pep_tag "ID=gene-" \
	--pep pep.fa --gff gene.gff \
	--out_gff ../Csat.gff --out_pep ../Csat.pep --out_list list.tsv

rm cds.list mrna.list alias.txt list.tsv gene.gff pep.fa

for i in Graimondii Esalsugineum Cpapaya Tcacao Crubella Ghirsutum Dstrictus Mperfoliatum Salba Cviolacea Alyrata Cgrandiflora Bstricta Camplexicaulis Alinifolium Rislandica Itinctoria Lannua Lsativum Tarvense Cmaritima Dsophioides Esyriacum Iamara Spinnata Evesicaria Chispanica Gmustelinum Gtomentosum Gbarbadense Gdarwinii Sparvula Ahalleri CsativaCAM116; do
	pigz -dc ${i}*gene.gff3.gz | awk '$3=="mRNA"&&$9~/longest=1/{print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' >gene.gff3
	pigz -dc ${i}*protein_primaryTranscriptOnly.fa.gz >protein.fa
	head -1 gene.gff3
	grep "^>" protein.fa | head -1
	perl gff_pep.pl \
		--abbr "${i}" --pep_tag "Name=" \
		--pep protein.fa --gff gene.gff3 \
		--out_gff ${i}.gff --out_pep ${i}.pep --out_list list.tsv
	rm gene.gff3 protein.fa list.tsv
done
