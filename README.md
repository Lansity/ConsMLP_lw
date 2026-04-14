# ConsMLP_lw

`ConsMLP_lw` 是一个轻量化的多层次超图划分工具（Multilevel Hypergraph Partitioner），面向 FPGA/EDA 场景做了实用优化，支持：

- 直接 `k` 路划分（`-mode direct`）
- 递归二分划分（`-mode recursive`）
- 三种约束模式：
  - 平衡约束（默认）
  - 类型约束（`-types`）
  - XML 绝对容量约束（`-xml`）

## 1. 功能概览

核心流程：

1. Coarsening：逐层收缩超图，降低问题规模
2. Initial Partitioning：在最粗层进行多次初始划分试探
3. Uncoarsening + Refinement：逐层投影回细图并做增量优化

当前实现特点：

- C++11 实现，单可执行文件
- 仅支持 `cluster` coarsen 和 `gfm` refine（接口保留，参数会做严格校验）
- 初始划分支持 `rand | ghg | ghg_opt | all`
- 递归模式支持任意 `k`（不要求 2 的幂）

## 2. 目录结构

```text
.
├── src/                    # 源码
├── include/                # 头文件
├── benchmarks/             # 示例 benchmark（.hgr/.insts）
├── compare_runs/           # 对比脚本与 XML 示例
├── tools/                  # 辅助脚本（如 XML 生成）
├── CMakeLists.txt
└── AGENTS.md
```

## 3. 编译

要求：

- CMake >= 3.10
- 支持 C++11 的编译器（g++/clang++）

Release：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Debug：

```bash
cmake -S . -B build_debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build_debug -j
```

生成可执行文件：

- `./build/ConsMLP_lw`
- `./build_debug/ConsMLP_lw`

## 4. 快速开始

查看帮助：

```bash
./build/ConsMLP_lw
```

基础二分示例：

```bash
./build/ConsMLP_lw benchmarks/case1.hgr -k 2 -mode direct -init rand
```

递归 4 分区示例：

```bash
./build/ConsMLP_lw benchmarks/case1.hgr -k 4 -mode recursive -init rand -passes 10
```

输出分区文件：

```bash
./build/ConsMLP_lw benchmarks/case1.hgr -k 4 -output case1.part.4
```

## 5. 参数说明

常用参数：

- `-k <num>`：分区数，默认 `2`
- `-mode <direct|recursive>`：划分模式，默认 `direct`
- `-init <rand|ghg|ghg_opt|all>`：初始划分策略，默认 `rand`
- `-imbalance <f>`：不平衡系数，默认 `0.05`
- `-passes <num>`：refine 最大迭代轮数，默认 `10`
- `-seed <num>`：随机种子，默认 `42`
- `-threshold <num>`：coarsening 阈值，默认 `100`
- `-coarsen cluster`：coarsen 算法（当前仅支持 `cluster`）
- `-refine gfm`：refine 算法（当前仅支持 `gfm`）
- `-output <file>`：输出分区结果文件

约束相关：

- `-types <file>`：节点类型文件（每行一个类型）
- `-xml <file>`：XML 容量约束文件（每个分区/SLR 的资源上限）
- `-relaxed <f>`：对 DSP/BRAM/IO 的放宽倍数，默认 `3.0`

可选开关：

- `-coarsen_opt`：启用 coarsen 过程中的性能优化策略

## 6. Types / XML 模式

### 6.1 Types 文件格式

每行一个节点类型，支持：

`LUT`, `FF`, `MUX`, `CARRY`, `IO`, `DSP`, `BRAM`, `OTHER`

示例：

```text
FF
LUT
DSP
```

### 6.2 XML 约束格式

按分区（SLR）定义资源上限，简化格式示例：

```xml
<SLR0>
<LUT><150040>
<FF><293483>
<DSP><1428>
<IO><40>
</SLR0>
<SLR1>
<LUT><150040>
<FF><293483>
<DSP><1428>
<IO><40>
</SLR1>
```

说明：

- XML 中未显式给出的类型默认视为无限制
- XML 包含非 LUT 资源约束时，必须同时提供 `-types`
- 在容量明显过剩时，程序会自动减少“激活分区数”，未激活分区保持空分区

### 6.3 自动生成 XML

项目提供脚本：

```bash
tools/gen_xml_from_types.sh <types_file> <output_xml> [--parts N] [--logic-util U] [--relaxed-extra R]
```

示例：

```bash
tools/gen_xml_from_types.sh \
  benchmarks/dsp_7400_lut_1110000_fd_2220000.insts \
  compare_runs/xml/dsp7400_4p_u80.xml \
  --parts 4 --logic-util 0.80 --relaxed-extra 0.15
```

## 7. 输出与日志

程序会在标准输出中给出：

- 超图统计（nodes/nets/pins、度分布）
- coarsening 各层规模
- 初始划分试验结果
- refine 改善情况
- 最终质量指标（Cut/Connectivity/SOED/Imbalance）

`-output` 文件中每行一个分区 ID，对应输入节点顺序。

## 8. 常见问题

1. `Unknown or incomplete option`
   参数拼写错误或缺少值，重新检查命令行。

2. `Invalid -mode / -init / -coarsen / -refine`
   当前版本只接受受支持的固定选项组合。

3. XML 模式报错要求 `-types`
   你的 XML 约束了非 LUT 资源，需要提供类型文件。

4. 推送到 GitHub 时认证失败
   请配置 SSH key 或 HTTPS 凭据（PAT）。

## 9. 许可与贡献

当前仓库未附带独立 License 文件。若用于团队协作，建议补充 `LICENSE` 与贡献规范（`CONTRIBUTING.md`）。
