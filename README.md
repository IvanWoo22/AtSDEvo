# AtSDEvo：拟南芥片段重复鉴定与 P/D 非对称序列演化

本项目用于鉴定拟南芥（*Arabidopsis thaliana*）片段重复
（segmental duplication, SD）事件，利用外群共线性区分祖先位置保留拷贝
P 与派生位置拷贝 D，并以复制前外群同源片段 P0 为祖先参考，
比较 P0→P 和 P0→D 两条分支上的派生序列变化。

## 当前主要结果

当前权威结果为 2026-07-29 完成的不依赖发生时间分箱的 P/D 分析：

- 1,655 个非冗余二拷贝 SD 事件；
- 754 个事件在 9 种区间与重叠阈值设置下保持 P/D 方向稳定；
- 239 个事件至少恢复 1 个事件匹配的 P0；
- 195 个事件进入三比对器 SNP 主终点，共 236,772 个可比位点；
- P 特异 SNP 9,498 个，D 特异 SNP 15,616 个，D/P = 1.644；
- 事件 bootstrap 的 D/P 95% 置信区间为 1.450–1.879；
- 1–10 bp 主证据级 micro-indel 为 P=23、D=54。

独立 RepeatMasker 全流程得到 237 个 SNP 主终点事件，P=11,723、D=19,912、
D/P=1.699；与 GFF 扩展主分析的 D/P=1.644 接近。跨掩膜 RO≥0.50 且均通过
age-free P/D 的 585 对事件中，P 位点方向一致率为 100%。Micro-indel 仍为
D>P，但效应量对掩膜和事件边界更敏感。

结果支持 D 分支积累了更多**派生序列变化**。这不等同于 D 具有更高的
原始突变率，因为选择、基因转换、DNA 修复、功能组成及 P0 可恢复性均可能影响
观察结果。

详见[当前结果报告](docs/results/current-age-free-results.md)和
[RepeatMasker 全流程重分析](docs/results/repeatmasker-full-reanalysis.md)，以及
[分析流程概览](docs/workflow-overview.md)。本地从冻结 MSA 重算可以精确恢复主
结果；重新执行比对器会出现少量平台相关边界差异，详见
[本地复现审计](docs/reproducibility-audit-2026-08-02.md)。

## 目录结构

```text
AtSDEvo/
├── README.md
├── environment.yml
├── docs/
│   ├── workflow-overview.md       # 整体流程
│   ├── reproducibility.md         # 复现命令
│   ├── project-context.md         # 科学背景与分析决策
│   ├── methodology/               # 方法细节
│   ├── results/                   # 冻结结果与版本比较
│   └── history/                   # 历史敏感性分析
├── workflow/                           # 按数据流编号的脚本
└── legacy/benchmark_v1/                # 早期 TAIR10/Col-CEN 原型
```

仓库不跟踪基因组大文件、BLAST 数据库、逐原子区比对缓存及完整结果树。
详见[复现说明](docs/reproducibility.md)和[数据管理原则](docs/data-policy.md)。

## 软件依赖

完整流程使用 RepeatMasker、BISER、BLAST+、MCScanX、MAFFT、MUSCLE5、PRANK，以及用于
增强型祖先状态分析的 IQ-TREE3。Python 汇总依赖 NumPy、SciPy 和 Matplotlib。
外部程序默认从 `PATH` 解析；也可用 `BISER_BIN`、`MCSCANX_BIN`、
`BLASTP_BIN`、`MAKEBLASTDB_BIN`、`GFFREAD_BIN` 和 `PYTHON_BIN` 显式指定。

仓库提供名为 `atsdevo` 的 Conda 环境：

```bash
conda env create -f environment.yml
conda activate atsdevo
python workflow/00_environment/check_environment.py --scope downstream
```

`environment.yml`安装 Python 分析包、BLAST+、MCScanX、三种序列比对器、
IQ-TREE、gffread、bedtools 与 samtools。HCC 的 BISER 1.4 Conda 包仅提供 Linux
x86-64 构建，并要求 Python 3.12，因此阶段 03 使用独立的
`environment-biser-linux.yml`，在 Linux 服务器或 `linux/amd64` 容器中运行。
已生成的 BISER calls 可直接用于 macOS 上的后续阶段。
RepeatMasker 只在重新生成独立 RM 掩膜时需要；已有冻结的 `.masked` 与 `.out`
文件时不影响主环境和下游复算。

## 术语

- **P**：保留在祖先共线位置的拷贝；
- **D**：位于派生基因组位置的拷贝；
- **P0**：复制前外群同源位置或祖先参考；
- **callable**：在物理坐标、序列大小写、gap、局部比对和祖先状态上都可判定的
  P/D/P0 位点。

P 不是比 D “更老的序列”；复制后的 P 和 D 分支经历了相同的时间。
