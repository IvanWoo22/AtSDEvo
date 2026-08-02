# RepeatMasker 掩膜敏感性分析

## 目的

片段重复识别依赖输入序列中哪些区域以小写形式被 BISER 忽略或降权。为判断
P/D 非对称结果是否由单一掩膜定义驱动，分析比较三种 TAIR12 soft-mask，并对
独立 RepeatMasker（RM）版本完整重算 BISER 及所有 P/D 下游步骤。

## 三种掩膜

| 名称 | 定义 | masked bp | 比例 |
|---|---|---:|---:|
| source | TAIR12 归档 FASTA 原有小写区域，生成方法不明 | 39,487,508 | 27.714% |
| GFF 扩展 | source 小写区域并入 GFF3 的 repeat/mobile-element 区间 | 51,953,668 | 36.464% |
| RepeatMasker | 对全大写 TAIR12 重新运行 RepeatMasker | 37,623,317 | 26.406% |

三份 FASTA 去除大小写后逐染色体完全一致。GFF 扩展与 RM 的交集为
31,304,755 bp，Jaccard 为 0.537；GFF-only 为 20,648,913 bp，RM-only 为
6,318,562 bp。因此两者不是“严格程度”不同的嵌套版本，而是来源和类别偏好不同
的掩膜定义。

RM 使用 RepeatMasker 4.2.3、RMBlast 2.14.1+ 和
`CONS-Dfam_withRBRM_3.9`，命令为：

```bash
RepeatMasker -xsmall -species arabidopsis -pa 6 TAIR12.uppercase.fna
```

独立 RM 新增区域以 rDNA/rRNA 为主；GFF 扩展对经典 TE 类别的覆盖更充分。
source 掩膜因来源未知，仅作为历史输入对照，不能解释为 RepeatMasker 或 DUST
结果。

## 比较层级

1. 验证三份 FASTA 的非大小写序列一致性，并比较 mask bp、交集和重复类别。
2. 在相同 BISER 1.4 默认阈值和 6 线程下独立识别 SD。
3. 用双臂最小 reciprocal overlap（RO）比较 calls；RO≥0.50 和 RO≥0.80
   分别表示中等和高边界一致性。
4. RM calls 独立执行年龄盲化网络规范化、外群共线性 P/D 极化、P0 映射、
   三比对器 SNP 和 de novo micro-indel 分析。
5. 比较完整端点，同时比较 RO≥0.80 的跨掩膜一对一事件，区分事件集合变化与
   同源事件内部的结果变化。

## 主分析选择

GFF 扩展继续作为主分析，因为它显式保留 TAIR12 原有掩膜，并补充参考注释中的
repeat/mobile-element；RM 全流程作为掩膜定义敏感性分析。两者应分别报告，不能
把 RM 新增事件直接并入主事件集。若后续评估混合方案，应预先定义为“RM 与从全
大写序列提取的 GFF-selected repeats 之并集”，避免再次继承来源未知的 source
小写区域。

## 对应脚本

- `workflow/01_reference/prepare_repeatmasker_reference.py`
- `workflow/01_reference/compare_three_softmasks.py`
- `workflow/03_sd_discovery/run_repeatmasker_biser.sh`
- `workflow/03_sd_discovery/compare_three_mask_biser_runs.py`
- `workflow/12_age_free_analysis/compare_mask_reanalyses.py`
