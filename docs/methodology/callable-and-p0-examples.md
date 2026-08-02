# Callable 与 event-matched P0 候选实例

## Callable

在当前直接 SNP endpoint 中，一个 callable site 是一个可以公平比较
`P0→P` 与 `P0→D` 的同源比对列。它必须同时满足：

1. P、D 和边界 P0 均有 A/C/G/T；
2. TAIR12 P、D 原始碱基均为大写；
3. 不位于长度≥11 bp的大gap及其两侧10 bp；
4. 局部半径25 bp内至少20个P/D成对碱基，局部错配≤40%；
5. MAFFT、MUSCLE5和固定树PRANK恢复相同P/D基因组坐标；
6. 可观察的更深P0不与边界P0冲突；
7. 重叠atom重复产生的同一P/D坐标只计一次。

例如 `SDJGI00057`：

- `Chr1:10156544(P)` ↔ `Chr1:12408881(D)`：P=A、D=A、P0=A，
  是 callable invariant；
- `Chr1:10156529(P)` ↔ `Chr1:12408866(D)`：P=A、D=T、P0=T，
  是 callable P-specific；
- `Chr1:10156599(P)` ↔ `Chr1:12408937(D)`：P=G、D=A、P0=G，
  是 callable D-specific。

该事件共有541个callable列，其中P-specific=19、D-specific=25。

## “无可接受 event-matched sequence candidate”

`event-matched` 表示事件X的P/D查询只能由事件X自己根据McScanX共线块生成的
外群候选区域接收。命中事件Y的候选区不会被归给事件X，以避免把其他旁系同源
或重复区域错误作为P0。

### 情况A：没有生成事件专属P0候选区

`SDJGI00107`，模式`3100`，P=locus_A：

- Alyrata状态3，成功生成P和D候选区；
- Bstricta状态1，是所需P0，McScanX记录了block126；
- 但block126中没有形成可用于该SD区段定位的锚点组合，因此没有
  `SDJGI00107|Bstricta|...`候选区；
- 最终记为`NO_SEQUENCE_CANDIDATE`。

这不表示Bstricta全基因组没有同源序列，只表示现行“先由共线锚生成候选区”
流程没有给该事件生成可搜索目标。

### 情况B：有候选区，但没有event-matched聚合HSP

`SDJGI00043`，模式`3202`，P=locus_B：

- Bstricta block43和Cviolacea block70均生成200 kb单侧锚候选区；
- 两个候选区的`two_sided_ordered_anchor_status=FAIL`；
- P查询没有在这些属于`SDJGI00043`的候选区中形成可汇总映射；
- 因而P-to-P candidate仍为NA。

### 情况C：有序锚失败

`SDJGI00016`，模式`3310`，P=locus_A：

- Dstrictus block27生成200 kb候选区；
- P序列identity=0.9565，但query coverage只有0.0331；
- 更重要的是候选区没有双侧有序锚，
  `two_sided_ordered_anchor_status=FAIL`；
- 即使存在一个高identity短命中，也不能证明它代表完整正交P0区段。

### 情况D：双侧锚通过，但覆盖不足

`SDJGI00046`，模式`3200`，P=locus_B：

- Bstricta block38生成57,271 bp候选区；
- 双侧有序锚通过；
- P查询identity=0.7490，但coverage=0.2051，低于0.35；
- 因而记为`QUERY_COVERAGE_FAIL`。

### 情况E：有全局同源命中，但目标属于其他事件

`SDJGI00310`，模式`3111`，P=locus_A：

- Bstricta、Dstrictus、Cviolacea均为状态1；
- 三个物种都没有为`SDJGI00310`生成事件专属P0候选区；
- 原始BLAST中可以看到其P/D命中Dstrictus序列，但目标FASTA记录属于
  `SDJGI00268|Dstrictus|...`、`SDJGI00211|Dstrictus|...`等其他事件；
- 现行汇总要求query event ID与target candidate event ID一致，因此这些命中
  被主动排除。

该实例直接说明`NO_SEQUENCE_CANDIDATE`不能解释为“外群无同源序列”，更准确的
含义是“没有在本事件预定义的共线候选区中形成可接受映射”。
