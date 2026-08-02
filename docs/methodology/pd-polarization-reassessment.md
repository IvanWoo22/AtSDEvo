# P/D 定向、过滤漏斗与 SNP 类型分析重审

日期：2026-07-29

## 1. 不以 time 为核心的 P/D 定向

四位数字依次表示 Alyrata、Bstricta、Dstrictus、Cviolacea：

- `0`：两个拟南芥位点均无 McScanX 共线支持；
- `1`：只有 locus_A 有共线支持；
- `2`：只有 locus_B 有共线支持；
- `3`：两个位点均有共线支持。

不要求单调 time 模式时，P/D 可直接按以下规则判定：

1. 至少出现一个 `1` 或 `2`；
2. 所有可用单拷贝状态必须支持同一物理位点，不能同时出现 `1` 和 `2`；
3. 在 representative/common-core/envelope × overlap 0.25/0.50/0.75
   的九种设置中，P 位点保持一致。

结果：

| 规则 | 事件 |
|---|---:|
| 1,655 个规范二位点事件 | 1,655 |
| 主模式可定向 P/D | 810 |
| 九种设置 P 稳定 | 754 |
| 其中至少两个单拷贝节点支持 | 344 |
| 其中只有一个单拷贝节点支持 | 410 |
| 原严格年龄/P-D 通过 | 585 |

754 个稳定 P/D 事件中，201 个是原 time 规则排除但可稳定定向的事件：
156 个原因为 `nonmonotonic_both_reappears`，45 个为
`boundary_node_uninformative`。time 可保留为描述性/敏感性变量，但无需再作为
第一部分 P/D 序列比较的主纳入条件。

建议把 344 个多节点支持事件设为 `PD-A`，410 个单节点但九设置稳定事件设为
`PD-B`；后者必须由实际 P0 序列映射进一步确认。

## 2. 过滤漏斗与基因组覆盖

“事件被过滤”与“事件没有观察到某类变异”必须分开。具有至少一个
micro-indel 的28个事件是结果产量，不是纳入过滤步骤，不应画在 eligibility
漏斗末端。

建议拆成：

1. **SD discovery/normalization**：BISER → network-stable calls →
   two-copy events；
2. **P/D eligibility**：age-free stable P/D → P0 序列可恢复 →
   continuous local atoms；
3. **assay QC**：三比对器一致、实际 callable 分母；
4. **biological outcomes**：SNP 数、micro-indel 数、无变异事件。

现有覆盖审计同时报告两拷贝区间长度之和及合并后的基因组 union：

| 阶段 | 数量 | 两拷贝区间和 | genome union |
|---|---:|---:|---:|
| BISER calls | 4,734 | 56.32 Mb | 21.54 Mb |
| network-stable calls | 1,865 | 11.72 Mb | 9.65 Mb |
| two-copy events | 1,655 | 10.30 Mb | 9.65 Mb |
| age-free stable P/D | 754 | 4.18 Mb | 4.09 Mb |
| 当前 time1–time3 | 362 | 2.06 Mb | 2.01 Mb |
| 当前至少一个 P0 | 110 | 0.681 Mb | 0.679 Mb |
| continuous atoms | 77 | 0.502 Mb | 0.502 Mb |
| SNP local MSA ≥200 | 70 | 0.470 Mb | 0.470 Mb |

第一步网络过滤删除60.6%的 calls，但 1,865→1,655 事件化只删除11.3%的记录，
同时几乎不损失 genome union（9.655→9.654 Mb），说明高比例冗余/复杂网络
过滤已经位于前端。

## 3. 362→110 的原因

当前“至少一个 P0”不是简单的 BLAST 命中条件。对于 time 指定的复制前
单拷贝物种，P 查询必须映射到该事件自己的 McScanX 候选区域，并通过：

1. event-matched candidate region；
2. 两侧有序共线锚点；
3. query coverage ≥0.35；
4. identity ≥0.55。

362 个事件中252个没有任何 P0 物种通过：

| 事件级失败概况 | 事件 |
|---|---:|
| 所有候选均无可接受 event-matched sequence candidate | 218 |
| 所有候选均失败双侧有序锚定 | 21 |
| 不同物种为混合失败原因 | 10 |
| 所有候选均 coverage 不足 | 3 |

252 个失败事件中，157 个在所有候选 P0 物种中根本没有 event-specific candidate
region；另95个至少有一个候选区域，但最终仍未通过映射。这里的
`NO_SEQUENCE_CANDIDATE` 不等于外群基因组中不存在同源序列，它可能只是：

- McScanX block 内没有位于 SD 两侧的可用基因锚；
- 候选区存在但 P 序列未形成 event-matched 聚合 HSP；
- 同源命中落到另一个事件的候选区，现行汇总会主动丢弃。

因此这一步是当前最值得重做的高损失环节。建议对754个 age-free P/D 事件重新
建立 P0 搜索：双侧锚候选为高置信层；单侧锚＋全基因组 P/D 联合定位为救援层；
要求 P 和 D 指向同一外群正交区域，而不是要求候选区必须先由双侧基因锚定义。

### `3111` 示例

`3111` 表示 Alyrata 两个位点都有共线支持，后三个物种只有 locus_A：
因此 locus_A=P、locus_B=D；time 解释是 time2，但 P/D 定向本身不依赖
“time2”标签。当前362事件中有15个 `3111`，6个获得 P0、9个失败；9个失败中
多数是三个 P0 物种均无 event-matched sequence candidate。

## 4. SBS96-like 三核苷酸谱

已按祖先 P0 碱基及其前后各1 nt 构建 96-channel spectrum；祖先为 A/G 时对
三核苷酸和替换共同反向互补，使中心统一为 C/T。只使用直接主 endpoint，
并按实际祖先三核苷酸 callable opportunity 归一化。

| 指标 | 数值 |
|---|---:|
| 直接 callable | 74,662 |
| 可恢复祖先三核苷酸 | 73,759 |
| P SNP | 2,033 |
| D SNP | 4,275 |
| P/D 原始谱 cosine similarity | 0.9759 |
| opportunity-normalized cosine similarity | 0.9748 |

初步结果显示 P、D 的三核苷酸谱形状高度相似，而总负担不同。当前更支持
“相似过程下 D 侧积累量增加”，尚不支持出现一个完全不同的新 mutation
signature。局部 channel 差异需要按事件 bootstrap/permutation 后再解释。

不建议直接把人癌 COSMIC signature 拟合结果解释为拟南芥机制。可采用其
SBS96 表示法，但机制参照应优先使用拟南芥 mutation-accumulation、甲基化、
UV/氧化损伤和 DNA repair 实验谱，并考虑 callable 区域的三核苷酸组成。

## 5. Transition/transversion 的正确拆分

互补归一后 transition 有两类：C>T、T>C；transversion 有四类：
C>A、C>G、T>A、T>G，不能把 transversion 再强行压成“两类”。

| 类型 | P | D | D/P |
|---|---:|---:|---:|
| C>T | 640 | 1,408 | 2.200 |
| T>C | 564 | 1,032 | 1.830 |
| C>A | 216 | 519 | 2.403 |
| C>G | 121 | 330 | 2.727 |
| T>A | 296 | 602 | 2.034 |
| T>G | 196 | 384 | 1.959 |

P 的 Ti/Tv=1.452，D 的 Ti/Tv=1.330；D/P 分别为 transition 2.027、
transversion 2.213。C>T/T>C 为 P=1.135、D=1.364。聚合位点检验提示
transition 内部组成可能不同，但碱基不是独立重复，正式结论应以事件为 block
进行重采样。

## 6. 功能区划分审查

当前脚本读取 TAIR12 GFF3 的所有 `gene`、`exon`、`CDS`，并按以下互斥优先级
逐个拷贝独立分类：

`CDS > exon_nonCDS > intron_or_gene_body > promoter_2kb > intergenic`

promoter 为每个 gene 按链方向上游2 kb；若与任何 gene 重叠，则 gene 类优先。
当前没有先冻结代表转录本，重叠基因/多转录本以“任一注释命中”形成 flags。

原表中 P intergenic callable=855、SNP=2，而 D intergenic
callable=2,297、SNP=85。虽然各自都使用了 callable 分母，但这两个集合不是
同一批同源位置：

- 855 个 P-intergenic 位点中，D 侧为 promoter 794、exon 30、
  intergenic 31；
- 真正 P/D 都为 intergenic 的只有31个位置，P=0、D=1；
- 只有1个事件同时具有 P/D intergenic 分母。

因此“P 的 intergenic 突变异常少”主要是功能背景不对称及样本极小，不能作为
生物学结论。主功能区分析应优先使用 P/D 同一 callable column 的
context-transition matrix，并把“双方同一 context”与“P/D 功能状态改变”
分开。

双方同一 context 时：

| context | callable pairs | P SNP | D SNP |
|---|---:|---:|---:|
| CDS | 46,591 | 1,359 | 2,936 |
| exon_nonCDS | 1,655 | 33 | 64 |
| intron/gene body | 7,270 | 292 | 536 |
| promoter | 511 | 11 | 15 |
| intergenic | 31 | 0 | 1 |

## 当前建议

1. time 降为协变量/敏感性分层，采用754个 age-free stable P/D 事件作为新的
   P0 搜索入口；
2. 以344个多节点支持事件为高置信层，410个单节点稳定事件为需序列验证层；
3. 重写 P0 候选生成，增加单侧锚和全基因组 P/D 联合定位救援；
4. 报告 eligibility funnel 时同时给事件数、两拷贝区间和、genome union；
5. micro-indel“有变异事件数”移出过滤漏斗，放到 outcome panel；
6. SBS96-like 谱先做事件 block-bootstrap，再与植物实验谱比较；
7. 功能区主结果改成 homologous-column paired context，冻结代表转录本规则。
