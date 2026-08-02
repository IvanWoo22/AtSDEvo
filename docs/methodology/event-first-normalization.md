# 事件优先的 SD 分类与 P/D 定向

事件分类按以下顺序进行：

```text
BISER calls
  → 物理拷贝位点
  → 稳定二位点事件
  → 外群共线状态
  → 发生时间区间与 P/D
```

早期 call-level 年龄标签只用于最终质控，不用于构建或纳入事件。

## 方法

1. 从 4,734 个 BISER calls 出发。
2. 在不使用任何年龄/P-D 证据的情况下建立 reciprocal overlap 0.5/0.8 网络，
   保留 1,865 个稳定 source calls 和 1,655 个二位点 components。
3. 在指定规范物理 `locus_A` 和 `locus_B` 之前，将每个 source-call arm 映射到网络 locus。
4. 选择代表 call，将两个 arm 重排到物理 A/B，同时计算 per-locus common core 与 envelope。
5. 在规范化事件区间上评估 4 个主外群物种的共线投影，然后才推断 `0/1/2/3`
   状态、年龄分箱和物理 P/D 位点。
6. 早期 1 kb 几何长度、uppercase、paired-M 和 joint-callable 指标只作质控，不决定纳入。
7. 用 overlap 0.25/0.50/0.75 和 representative/common-core/envelope 进行敏感性检查。

## 事件数量

- 规范化二位点事件：1,655；
- 不使用长度前置门槛时，事件优先年龄/P-D pass：585；
- time1=187、time2=71、time3=129、time4=198；
- 在 overlap 0.25/0.50/0.75 下分类一致：542；
- 阈值稳定的 time1–time3 候选：362。

362 个候选均含有至少 1 个 paired 且 jointly uppercase-callable 碱基，因此全部进入包容性
局部序列分析。其中 130 个事件的 jointly callable 长度低于 1 kb，但仍保留并使用
实际 callable 分母。

相对早期基线，不设长度门槛新增 139 个事件，没有事件仅因不满足 1 kb 而被删除。

## 主要输出

- `event_first_catalog.tsv`：不依赖年龄的物理事件目录；
- `event_first_events.tsv`：1,655 个事件级分类；
- `source_arm_to_event_locus.tsv`：call arm 到物理 locus 的显式映射；
- `event_first_threshold_stable.tsv`：阈值稳定年龄/P-D 事件；
- `event_first_time1_time3_threshold_stable.tsv`：time1–time3 候选集；
- `event_classification_by_scope_threshold.tsv`：范围与重叠阈值敏感性；
- `event_effective_unmasked_length_qc.tsv`：几何、uppercase、paired-M 和 joint-callable 长度审计。

## Time4/time5 的缺失容忍性枚举

另设一个系统发育约束的敏感性分类，允许单拷贝状态 `1/2` 代表谱系特异的复制后丢失。
在 ≥1 kb 条件下，扩展集包含 202 个 time4 事件（其中 183 个阈值/范围稳定）和
375 个 time5 事件（其中 365 个稳定）。Time5 比最远主节点更早，仅凭这 4 个物种无法得到
有限更老边界或可靠 P/D 极化，因此不进入当前主序列终点。
