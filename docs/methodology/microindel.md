# time1–time3 de novo micro-indel分析

## 分析定义

本分析不读取BISER CIGAR的`I/D`。BISER仅用于确定已经接受的SD P/D位点和
祖先同源核心坐标。连续P/D、近缘复制后P/D和复制前P0单倍型分别经过：

1. MAFFT L-INS-i de novo gap发现；
2. MUSCLE 5独立发现；
3. 固定物种树PRANK v.250331 `+F`精确坐标验证；
4. 重叠窗口坐标规范化和相邻同方向gap合并；
5. 外群状态极化、重复/soft-mask/低复杂度过滤。

time2/time3要求至少一个复制后物种复现P/D存在-缺失模式，并由P0判向。
time1没有复制后物种P/D，主终点改为边界P0与至少一个更深P0状态一致，且
所有可用P0无冲突。

## 筛选结果

| 阶段 | 数量 |
|---|---:|
| 具备所需外群状态的连续原子区域 | 1,245 |
| 具有MAFFT/MUSCLE共同gap的原子区域 | 292 |
| 三算法精确坐标一致、去重后的候选 | 78 |
| 可极化候选（含序列敏感性） | 44 |
| 高质量主调用 | 38 |
| 含主调用的独立SD事件 | 21 |

## 年龄分层

| 年龄 | SD事件 | P micro-indel | D micro-indel | D>P事件 | P>D事件 | 相等 | sign test p |
|---|---:|---:|---:|---:|---:|---:|---:|
| time1 | 31 | 7 | 19 | 11 | 2 | 18 | 0.02246 |
| time2 | 7 | 0 | 0 | 0 | 0 | 7 | NA |
| time3 | 17 | 1 | 11 | 7 | 0 | 10 | 0.01563 |
| 合计 | 55 | 8 | 30 | 18 | 2 | 35 | 0.000402 |

总体D/P调用数为30/8=3.75。调用数检验仅作描述；主推断采用SD事件作为
独立单位。

具体方向构成为：

- time1：P insertion 3、P deletion 4、D insertion 6、D deletion 13；
- time3：P insertion 1、P deletion 0、D insertion 1、D deletion 10；
- 合计：P insertion 4、P deletion 4、D insertion 7、D deletion 23。

因此D侧增加主要来自deletion，但insertion也保持D>P。该拆分为描述性结果，
不能把单个indel调用当作相互独立重复。

time1的26个主调用证明多P0方案能够补充原先缺失的年龄层。time2仍无主调用，
原因是只有一个复制后验证物种且候选复现不足；不能把0解释为time2没有indel。

## 窗口敏感性

替代参数使用80 bp junction锚点、300 bp窗口和250 bp步长：

- 主参数38个主调用，替代参数37个；
- 28/38个主调用精确坐标复现；
- 14个事件在两套参数中均D>P，2个均P>D；
- 稳健事件方向sign test：p=0.00418。

分层稳健方向为：

- time1：10个D>P、2个P>D，p=0.03857；
- time2：无可检验事件；
- time3：4个D>P、0个P>D，p=0.125。

因此总体D偏向具有跨窗口参数支持；time1单独达到提示性/统计支持，time3方向
一致但稳健事件数偏少。

## 结论边界

当前结果支持“D分支积累更多短indel”的先决证据，并与SNP的D偏高方向一致。
但应称为派生变化积累差异，而不是绝对mutation rate差异；选择、基因转换、
局部组装质量和外群可比性仍可能共同影响结果。

## 主要文件

- `denovo_microindel_inference.tsv`
- `age_stratified_PD_microindel_summary.tsv`
- `event_level_PD_microindel_rates.tsv`
- `event_level_paired_statistical_summary.tsv`
- `window_parameter_call_robustness.tsv`
- `window_parameter_event_robustness.tsv`
- `window_parameter_age_robustness.tsv`
- `alignments/PRANK_FIXED_TREE/`（含祖先和事件输出）
