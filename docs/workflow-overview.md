# 当前分析流程概览

## 数据流

```text
TAIR12 assembly + GFF3
  → 染色体命名标准化与 soft-mask 审计
  → 注释扩展 soft-mask
  → BISER：4,734 个局部 calls
  → 年龄盲化的 call/locus 网络规范化
  → 1,655 个非冗余二拷贝事件
  → 4 个外群物种的 MCScanX 状态
  → 9 种设置下 age-free P/D 稳定：754 事件
  → 事件匹配的外群 P0 映射：239 事件
  → 连续 P/D/P0 原子区：208 事件
  ├─ 局部 SNP：MAFFT + MUSCLE5 + PRANK → 195 事件主终点
  └─ de novo 1–10 bp indel：三比对器独立 gap 发现 → 77 个主调用
```

## 定义当前终点的关键决策

1. BISER call 不直接当作独立生物学事件，须先规范化为稳定二拷贝网络。
2. P/D 方向不依赖单调的复制时间分箱。至少 1 个单拷贝状态必须支持 P，
   所有可观察单拷贝状态必须一致，且 P 在 3 种事件范围和 3 个重叠阈值下不变。
3. MCScanX 状态 `0` 表示未检出，不表示已证实的生物学缺失。
4. 早期 1 kb 整事件条件仅作诊断列，不作为前置纳入门槛。
5. P0 必须与当前事件匹配，并通过有序锚点、query coverage 和 identity 阈值。
6. SNP 只保留 MAFFT、MUSCLE5 和固定树 PRANK 在 P/D 物理坐标、祖先碱基和
   分支类别上完全一致的位点。
7. Micro-indel 由局部单倍型重新比对后 de novo 发现，不把 BISER CIGAR I/D 当作变异调用。
8. 统计推断单位是 SD 事件，而非单个变异位点。

## 主终点与敏感性分析

- SNP 主终点：age-free、三比对器一致、每事件至少 200 个 callable sites。
- Micro-indel 主终点：age-free P0 队列中的三比对器高质量 1–10 bp calls。
- 保守验证：边界 P0 与更深 P0 直接一致。
- 增强型敏感性：固定拓扑 JC/K2P 祖先状态重建。
- 历史对照：44/51/55/70/71 事件的严格投影或 time-gated 终点。

详细定义见 `methodology/`，结果变化见 `results/result-evolution.md`。

