# TAIR12 片段重复溯源与 P/D 非对称演化

最后更新：2026-07-29

## 1. 科学问题

本项目研究完整片段重复（segmental duplication, SD），而非仅研究重复基因。
当前分析聚焦于第一个先决问题：能否在 TAIR12/Col-0 中建立非冗余 SD 事件，
稳定区分 P 和 D，并比较共同祖先状态 P0 到 P、D 两条分支的派生序列变化。

后续研究可在此基础上扩展至：

1. 多个拟南芥生态型中的 PAV/CNV、结构状态、种内多样性和基因转换；
2. Col-0 多组学中 P/D 的 5mC、组蛋白修饰、染色质可及性、小 RNA 和表达差异；
3. copy-resolved 甲基化数据的独立验证；
4. 跨时间尺度比较突变、选择、基因转换和表观状态。

## 2. 核心术语与解释边界

- **P**：保留在祖先共线位置的拷贝；
- **D**：位于派生基因组位置的拷贝；
- **P0**：复制前外群同源位置或祖先参考；
- **branch-specific derived change**：相对 P0 可极化到 P 或 D 分支的派生变化。

P 不是“更老的序列”，D 也不是“更年轻的序列”。两条复制后分支经历了相同时间。
当前测量反映派生变化的积累，不直接等同于原始突变率。

## 3. 参考组与 soft-mask

- 参考组为 TAIR12/Col-0；
- 染色体名称统一为 `Chr1`–`Chr5`；
- 比较源 FASTA 已有 soft-mask 与注释扩展 soft-mask；
- 正式 BISER 输入使用注释扩展 soft-mask；
- 源小写区域不能等同于当前 GFF 重复注释的完整集合。

## 4. SD 识别与事件化

BISER 1.4 在正式 mask 上产生 4,734 个局部 calls。BISER call 是局部比对记录，
不是非冗余生物学事件，因此分析采用：

```text
BISER calls
  → 年龄盲化 call/locus 网络
  → reciprocal overlap 0.5 和 0.8 下拓扑稳定
  → 1,865 个稳定 source calls
  → 1,655 个非冗余二拷贝事件
```

每个事件保留 representative interval、per-locus common core 和 envelope，用于范围敏感性分析。

## 5. 外群共线性

使用 4 个主节点物种：

| 缩写 | 物种 | 数据版本 |
|---|---|---|
| Alyrata | *Arabidopsis lyrata* | v2.1 |
| Bstricta | *Boechera stricta* | v1.2 |
| Dstrictus | *Diptychocarpus strictus* | v2.1 |
| Cviolacea | *Cleome violacea* | v2.1 |

成对 BLASTP + MCScanX 统一使用 `-s 3 -m 2 -w 0`。每个物种对一个二拷贝事件的状态为：

- `0`：两个位点均未检出达阈值的共线支持；
- `1`：仅 locus_A 有支持；
- `2`：仅 locus_B 有支持；
- `3`：两个位点均有支持。

`0` 代表未检出，不代表已证实缺失。

## 6. Age-free P/D 定向

发生时间分箱对研究 SD 年龄有用，但不再决定事件能否进入 P/D 序列比较。
正式 P/D 条件为：

1. 4 物种状态中至少出现 1 个 `1` 或 `2`；
2. 所有可观察的 `1/2` 状态指向同一个拟南芥位点；
3. P 位点在 representative/common-core/envelope 三种范围与
   0.25/0.50/0.75 三个重叠阈值下完全一致。

共得到 754 个稳定 P/D 事件：344 个有至少 2 个单拷贝节点支持，410 个仅有
1 个单拷贝节点支持，但在 9 种设置下方向稳定。

## 7. P0 恢复

所有与 P 方向一致的单拷贝物种均可作为 P0 候选。映射要求：

- query event ID 与 target candidate event ID 一致；
- 通过两侧有序共线锚点；
- query coverage ≥ 0.35；
- 非重叠加权 identity ≥ 0.55。

BLAST HSP 召回使用 `evalue=1e-3`。E-value 只用于候选召回，最终接受仍由固定
coverage、identity 和锚点规则决定。754 个事件中 239 个至少恢复 1 个 P0，其中
167 个只有 1 个 P0 物种，72 个有至少 2 个 P0 物种。

## 8. SNP 方法

连续 P/D/P0 片段被拆分为不跨越映射断点的原子区，分别运行 MAFFT L-INS-i、
MUSCLE5 和固定树 PRANK。位点须同时满足：

- 三种比对恢复相同 P/D 物理坐标、祖先碱基和分支类别；
- 去除 P–D、P–P0 或 D–P0 的 >10 bp gap 及两侧 10 bp；
- 局部±25列至少有20个 P/D 配对碱基，局部错配率≤40%；
- P、D、P0 都是非 soft-mask A/C/G/T；
- 更深 P0 若可观察，不得与边界 P0 冲突。

主终点要求每事件至少 200 个三比对器一致位点。

## 9. Micro-indel 方法

正式 1–10 bp micro-indel 不读取 BISER CIGAR I/D 作为调用。BISER 只提供 SD 位点和
同源坐标，局部 P/D/P0 单倍型由 MAFFT、MUSCLE5 和固定树 PRANK 重新比对，
候选 gap 必须在三种算法中坐标一致，并通过大小写、注释重复、低复杂度和 P0 方向质控。

## 10. 当前结果

### SNP

- 195 个事件；
- 236,772 个 callable sites；
- P=9,498，D=15,616，D/P=1.644；
- 148 个事件 D>P，43 个 P>D，4 个相等；
- sign test `p=1.07e-14`；
- paired Wilcoxon `p=1.63e-18`；
- 事件 bootstrap D/P 95% CI：1.450–1.879。

### Micro-indel

- 239 个 P0 队列事件；
- 77 个主调用：P=23，D=54，D/P=2.348；
- 44 个事件有主调用；
- 事件 sign test `p=0.00643`；
- paired Wilcoxon `p=0.02797`；
- bootstrap D/P 95% CI：1.167–5.300。

Micro-indel 方向与 SNP 一致，但调用更稀疏、置信区间更宽，因此证据等级较低。

## 11. 主要限制

- P0 队列受外群注释、共线性和事件锚定可恢复性影响；
- 两侧锚点提高事件特异性，但可能损失基因间或重排事件；
- 410 个 single-node P/D 事件的证据弱于 344 个 multi-node 事件；
- 三比对器一致性可降低 alignment artifact，但不是祖先状态的绝对真值；
- P/D 功能背景不是随机分配，必须优先比较同源位点的配对功能背景；
- D>P 可由突变输入、选择、修复、基因转换、局部重组或保留偏倚共同造成。

## 12. 优先后续分析

1. 比较 344 个 multi-node 事件与 410 个 single-node 事件的 SNP 效应；
2. 建模评估 754→239 的 P0 可恢复性偏倚；
3. 对事件长度、repeat 比例和功能背景进行匹配或加权；
4. 系统筛查基因转换；
5. 进行 P0 物种 leave-one-out 与更深 P0 冲突敏感性；
6. 在 239 事件集上补充 micro-indel 替代窗口参数复算。

