# 分析流程与脚本分工

编号目录按数据流排列，而非按脚本编写时间排列。

| 阶段 | 职责 | 主要产物 |
|---|---|---|
| 01 | TAIR12 参考组与注释扩展 soft-mask | 标准化 FASTA/GFF、mask 区间 |
| 02 | 源 mask 与注释扩展 mask 覆盖审计 | 覆盖统计 |
| 03 | BISER 识别与 mask 策略比较 | SD calls、资源使用日志 |
| 04 | 4 个主外群物种的资源准备 | genome/GFF/peptide 输入 |
| 05 | 成对 BLASTP 与 MCScanX | 外群共线区块 |
| 06 | call-to-event 规范化与 P/D 状态审计 | 1,655 个二拷贝事件 |
| 07 | 事件匹配 P0 候选生成与映射 | P0 队列与映射结果 |
| 08 | 纳入条件与高错配敏感性 | 历史对照终点 |
| 09 | 局部 MSA SNP 与 de novo micro-indel | 事件级和位点级指标 |
| 10 | 图表与 PDF 输出 | 报告素材 |
| 11 | 固定拓扑祖先状态敏感性 | JC/K2P ASR 汇总 |
| 12 | age-free P/D 重评与最终汇总 | 当前 195 事件主终点 |

脚本在独立分析目录中读写大型输入、中间结果和比对缓存。执行顺序见
[`docs/reproducibility.md`](../docs/reproducibility.md)。阶段 01–06 是上游前提；
阶段 12 是当前主分析；阶段 08 和 11 为敏感性分析。

