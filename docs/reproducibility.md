# 复现说明

代码仓库与大型分析数据目录分离。以下约定：

```bash
CODE=/path/to/AtSDEvo
P=/path/to/TAIR12_outgroup_comparative_genomics
export SD_PROJECT_ROOT="$P"
export SD_AGE_FREE_ROOT="$P/15_age_free_pd_sequence_variation"
```

`P`中存放参考组、外群资源、中间结果和比对缓存；`CODE`中存放本仓库追踪的
脚本。完整上游数据须先按`workflow/01`至`workflow/07`准备。

## 1. 不依赖时间分箱的事件与 BISER 同源 core

```bash
python3 $CODE/workflow/12_age_free_analysis/build_age_free_event_input.py
python3 $CODE/workflow/07_p0_mapping/build_high_priority_cigar_cores.py \
  --project $P \
  --output $P/15_age_free_pd_sequence_variation \
  --events $P/15_age_free_pd_sequence_variation/inputs/age_free_pd_events.tsv \
  --age-free --minimum-core-bp 1 --minimum-callable-bp 1
python3 $CODE/workflow/07_p0_mapping/build_outgroup_candidate_regions.py \
  --project $P --pilot $P/15_age_free_pd_sequence_variation
```

## 2. BLAST 候选映射

每个`SPECIES`为`Alyrata Bstricta Dstrictus Cviolacea`。先对各候选FASTA运行
`makeblastdb`，再运行：

```bash
blastn -task blastn \
  -query $P/15_age_free_pd_sequence_variation/core/TAIR12_PD_homologous_cores.fa \
  -db $P/15_age_free_pd_sequence_variation/outgroup_mapping/blastn/SPECIES \
  -evalue 1e-3 -word_size 11 -dust no -soft_masking false \
  -max_target_seqs 100000 \
  -outfmt '6 qseqid sseqid qlen length pident mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out $P/15_age_free_pd_sequence_variation/outgroup_mapping/blastn/SPECIES.hits.tsv

blastn -task blastn \
  -query $P/15_age_free_pd_sequence_variation/core/TAIR12_PD_homologous_cores.fa \
  -db $P/15_age_free_pd_sequence_variation/outgroup_mapping/blastn/SPECIES \
  -evalue 1e-3 -word_size 11 -dust yes -soft_masking true \
  -max_target_seqs 100000 \
  -outfmt '6 qseqid sseqid qlen length pident mismatch gapopen qstart qend sstart send evalue bitscore qstrand sstrand qseq sseq' \
  -out $P/15_age_free_pd_sequence_variation/outgroup_mapping/blastn_aligned/SPECIES.aligned_hits.tsv
```

本机BLAST 2.17.0不输出不支持的`qstrand` token，实际aligned结果为16列，
解析器把query设为plus，field 14作为subject strand。

```bash
python3 $CODE/workflow/07_p0_mapping/summarize_candidate_blastn.py \
  --pilot $P/15_age_free_pd_sequence_variation \
  --min-query-coverage 0.35 --min-identity 0.55
python3 $CODE/workflow/12_age_free_analysis/build_age_free_p0_queue.py
```

## 3. Micro-indel 与原子区

```bash
python3 $CODE/workflow/09_sequence_variation/analyze_denovo_msa_microindels.py \
  --project $P \
  --pilot $P/15_age_free_pd_sequence_variation \
  --inclusive-queue $P/15_age_free_pd_sequence_variation/pilot/age_free_p0_mapping_queue.tsv \
  --output $P/15_age_free_pd_sequence_variation/microindel_local_msa \
  --age-bins age_free --event-set mapping_union --age-free-pd
```

## 4. SNP 局部 MSA 与类型

```bash
python3 $CODE/workflow/09_sequence_variation/analyze_local_msa_snps.py \
  --project $P \
  --pilot $P/15_age_free_pd_sequence_variation \
  --inclusive-queue $P/15_age_free_pd_sequence_variation/pilot/age_free_p0_mapping_queue.tsv \
  --atom-root $P/15_age_free_pd_sequence_variation/microindel_local_msa \
  --output $P/15_age_free_pd_sequence_variation/snp_local_msa \
  --prank-cache $P/15_age_free_pd_sequence_variation/microindel_local_msa/alignments/PRANK_FIXED_TREE \
  --p0-rule boundary_no_conflict --prank-event-threshold 200

python3 $CODE/workflow/09_sequence_variation/analyze_local_msa_snp_types.py \
  --project $P \
  --sites $P/15_age_free_pd_sequence_variation/snp_local_msa/three_aligner_local_callable_sites.tsv \
  --events $P/15_age_free_pd_sequence_variation/snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv \
  --output $P/15_age_free_pd_sequence_variation/snp_types_ge200

python3 $CODE/workflow/12_age_free_analysis/analyze_snp_context_and_function.py \
  --analysis-root $P/15_age_free_pd_sequence_variation \
  --output $P/15_age_free_pd_sequence_variation --minimum-sites 200

python3 $CODE/workflow/12_age_free_analysis/summarize_age_free_results.py
```

## 5. 软件

- BISER 1.4：原运行位于`hvnlr`环境；
- MCScanX：`mcscan`环境，运行参数`-s 3 -m 2`；
- NCBI BLAST+ 2.17.0；
- MAFFT 7.526；
- MUSCLE 5.3；
- PRANK 250331；
- IQ-TREE 3.1.1（历史增强ASR使用，当前直接P0主结果不依赖）；
- Python分析环境需NumPy、SciPy、Matplotlib；
- PDF：Pandoc 3.10.1、XeLaTeX/TeX Live 2024。

仓库中的当前脚本已移除旧服务器绝对可执行路径。BISER、MCScanX/BLAST、
PRANK和IQ-TREE默认从`PATH`解析，也可通过脚本参数或`*_BIN`环境变量显式指定。
