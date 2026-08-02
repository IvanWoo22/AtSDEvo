# RepeatMasker 全流程重分析

完成日期：2026-08-03

## 分析设置

- BISER：1.4，默认阈值，6 线程；
- 事件：年龄盲化网络、最小长度 1 kb、RO 阈值 0.25/0.50/0.75；
- P/D：representative、common core、envelope 三种范围均指向相同 P；
- P0：BLASTN `evalue=1e-3`、`word_size=11`、query coverage ≥0.35、
  identity ≥0.55；
- SNP：MAFFT L-INS-i、MUSCLE 5.3、PRANK 250331 三算法一致，
  `boundary_no_conflict`，每事件至少 200 callable sites；
- micro-indel：age-free、mapping-union、三比对器一致。

并行预计算只改变任务调度，不改变比对参数。4,562 个原子区的 PRANK 任务中
4,557 个成功，5 个失败原子区被排除。

## 全流程漏斗

| 阶段 | GFF 扩展 | RepeatMasker | RM/GFF |
|---|---:|---:|---:|
| BISER calls | 4,734 | 6,669 | 140.87% |
| 严格二拷贝事件 | 1,655 | 1,892 | 114.32% |
| 稳定 age-free P/D | 754 | 864 | 114.59% |
| 至少一个映射 P0 | 239 | 262 | 109.62% |
| 原子区 | 4,457 | 4,562 | 102.36% |
| SNP ≥200 事件 | 195 | 237 | 121.54% |

## 端点比较

| 端点 | GFF 扩展 | RepeatMasker |
|---|---:|---:|
| SNP P | 9,498 | 11,723 |
| SNP D | 15,616 | 19,912 |
| SNP D/P | 1.644 | 1.699 |
| micro-indel P | 23 | 20 |
| micro-indel D | 54 | 58 |
| micro-indel D/P | 2.348 | 2.900 |

RM 版本的 SNP 三核苷酸谱与 GFF 版本高度相似：原始计数余弦相似度 0.9941，
按机会归一化后为 0.9928。结果说明掩膜主要改变事件发现和边界，而不是把 SNP
替换谱改造成另一种模式。

## 解释

完整端点同时包含事件集合、边界、P0 可恢复性与 callable 序列的变化；跨掩膜
一对一事件更接近同一生物学对象的直接比较。两种比较均支持 SNP 的 D>P，因而
该结论具有较强掩膜稳健性。micro-indel 数量稀疏且更依赖 gap 边界，应报告方向
和敏感性范围，不把单一 D/P 数值作为固定效应量。
