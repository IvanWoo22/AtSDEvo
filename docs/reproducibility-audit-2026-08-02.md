# 2026-08-02 本地复现审计

## 测试环境

- 平台：macOS Apple Silicon（arm64）；
- Python 3.11.15；
- NumPy 2.4.6、SciPy 1.17.1、Matplotlib 3.11.1；
- MAFFT 7.526、MUSCLE 5.3、PRANK 250331、IQ-TREE 3.1.3；
- BLAST+ 2.16.0；
- gffread 0.12.9、bedtools 2.31.1、samtools 1.24；
- MCScanX 1.0.0：Bioconda macOS arm64 原生包；
- 测试数据：`DRMEGD/TAIR12_outgroup_comparative_genomics` 的写时复制副本。

所有测试输出均写入临时数据副本，未覆盖权威结果目录。
环境检查、BLAST 核酸库构建与短序列查询均实际运行通过。

## 可精确复现的结果

以下输出与同步结果逐字节一致：

- 年龄无关 P/D 事件输入：754 个事件；
- P0 映射队列：239 个事件；
- SNP 功能与三核苷酸背景全部 TSV；
- 使用冻结 MSA 文件重新计算的 SNP 汇总：195 个事件、236,772 个 callable
  位点、P=9,498、D=15,616；
- 使用冻结 MSA 文件重新计算的六类 SNP 谱和功能背景全部 TSV；
- 最终统计汇总中的所有 TSV。

图形文件的像素或 PDF 元数据可因字体缓存与渲染后端不同而变化，不作为数值复现
判据。

## 重新执行比对器时的差异

从原子序列开始重新运行 MAFFT、MUSCLE5 和 PRANK 后，方向性结论不变，但少量
边界位点发生变化：

| 终点 | 冻结结果 | 本机重新比对 |
|---|---:|---:|
| SNP 事件数 | 195 | 195 |
| SNP callable 位点 | 236,772 | 236,784 |
| P 特异 SNP | 9,498 | 9,513 |
| D 特异 SNP | 15,616 | 15,618 |
| SNP D/P | 1.644 | 1.642 |
| 主证据 micro-indel | P=23，D=54 | P=23，D=56 |
| 有主证据 micro-indel 的事件 | 44 | 46 |

差异位于多序列比对边界，而不是 P/D 事件定义、P0 队列或下游统计代码。相同主
版本的比对器在不同平台构建、浮点实现或启发式路径下仍可能产生少量不同 gap
布局。因此：

1. 当前论文级数值继续以冻结结果为准；
2. 复核下游统计时应复用冻结的原子 MSA；
3. 若从原始序列完全重跑，应把比对器平台和构建版本视为分析版本的一部分，并将
   新结果作为独立版本比较，不与冻结结果混合；
4. 不应把两条新增 D micro-indel 直接并入当前主结果，除非完成逐位点人工复核和
   预先定义的版本更新。

## 本次代码修正

- 将 Python 固定为 3.11，避免环境自动升级到未经测试的 3.14；
- 删除主流程未使用的 Biopython 依赖；
- 修正 PRANK 在 `PATH` 中存在但被误判为缺失的问题；
- 新增分范围环境检查，能够识别 MUSCLE 3.x 误用；
- 将 MCScanX 改为 Bioconda 原生依赖；
- 为 HCC 的 Linux x86-64 BISER 1.4 增加独立环境文件。

## Conda 包复核

- Bioconda `mcscanx=1.0.0` 在独立 macOS arm64 环境中安装成功；可执行文件为
  Mach-O arm64，并以最小 GFF 和 BLAST 输入完成真实运行、生成 collinearity 与
  HTML 输出。
- HCC `biser=1.4` 在 macOS arm64 求解时返回 `PackagesNotFoundInChannelsError`。
  该包的唯一构建为 Linux x86-64 ELF，强制下载后在 macOS 执行返回状态 126
  (`exec format error`)。因此它可用于 Linux 服务器或 `linux/amd64` 容器，不能
  作为 macOS 原生依赖。`environment-biser-linux.yml` 已按 `linux-64` 完成 Conda
  求解校验，能够解析到 HCC `biser=1.4=py312h06f12e4_0` 及其全部依赖。
