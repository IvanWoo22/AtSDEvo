# 4 个 JGI 主节点物种的 MCScanX 共线性分析

## 输入与运行有效性

所有 JGI 文件均通过 manifest MD5 校验。*A. lyrata*、*B. stricta*、*D. strictus* 和
*C. violacea* 分别包含 31,073、27,416、27,066 和 21,850 条代表蛋白。
每个 peptide ID 都对应 1 个 primary mRNA，所有 GFF 坐标均对应组装序列且未越界。
4 组 BLASTP 与 MCScanX 任务退出码均为 0，collinearity 输出中未发现未知 ID。

## 成对结果

| 主节点 | 共线 blocks | 基因对 | 唯一 TAIR12 基因 | TAIR12 覆盖率 |
|---|---:|---:|---:|---:|
| *A. lyrata* | 615 | 22,612 | 21,824 | 81.23% |
| *B. stricta* | 645 | 20,649 | 19,959 | 74.29% |
| *D. strictus* | 777 | 16,518 | 16,037 | 59.69% |
| *C. violacea* | 885 | 11,668 | 11,637 | 43.31% |

唯一 TAIR12 基因覆盖率随系统发育距离单调下降，符合预期节点顺序。原始 block 数向深层比较
反而增加，因此不能将 block 数当作时钟：碎片化或更短 block 可以在唯一锚点覆盖率下降时
使 block 数增加。

在相同 `-w 0` 参数下，*A. lyrata*、*D. strictus* 和 *C. violacea* 的汇总结果与
早期运行一致；*B. stricta* 新增 1 个基因对和 1 个唯一 TAIR12 基因。

这些共线 blocks 只用于判定 TAIR12 SD 两个拷贝是否与 TAIR12-projected 共线区域重叠；
不在外群中重新进行 SD discovery。

