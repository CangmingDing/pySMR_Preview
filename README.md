# pySMR_Preview

<div align="center">

**基于 Python 的SMR分析工具（v0.3.0 - 功能完整版）**

[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)

[English Version](README_EN.md) | 中文文档

</div>

---

## 📖 项目简介

**pySMR_Preview** 是对杨剑教授团队开发的原始 [SMR 软件](https://yanglab.westlake.edu.cn/software/smr/#Overview) 的 Python 实现（v0.3.0）。本项目严格遵循开源政策，**已完整复刻 SMR 和 HEIDI 测试的核心功能**，并在数据处理层面实现了更鲁棒的等位基因匹配策略。

> **✅ 功能对齐状态**  
> 当前版本已实现与 C++ 原版在**分析层面的完全对齐**：
> - ✅ 单 SNP SMR 测试（Top-SNP SMR）
> - ✅ **多 SNP SMR 测试（Multi-SNP / Set-based SMR）** - **新增功能！**
> - ✅ HEIDI 异质性检验
> - ✅ 更鲁棒的等位基因匹配（支持链翻转和互补碱基匹配）

本工具旨在为整合 GWAS 和 eQTL 摘要数据提供灵活、高效的 Python 解决方案。

---

## ✨ 核心功能

- ✅ **单 SNP SMR 测试**：使用单个遗传变异作为工具变量，检测基因表达与复杂性状之间的关联
- ✅ **🎯 多 SNP SMR 测试** *(新增)*：  
  整合区域内多个独立 eQTL 信号进行联合检验（Set-based test），使用 Satterthwaite 近似计算加权卡方和的 P 值，显著提升统计功效
- ✅ **HEIDI 测试**：检测依赖工具变量中的异质性，以区分多效性（Pleiotropy）和连锁（Linkage）
- ✅ **高级等位基因匹配**：自动处理链翻转和碱基互补（A/G ↔ T/C），最大化数据利用率
- ✅ **完整可视化**：生成 Locus 轨迹图、Effect Size 散点图和 **Manhattan 图**
- ✅ **发表级绘图**：支持 PDF 矢量输出，Helvetica 字体，以及 **Gene Symbol 自动转换**
- ✅ **高效 LD 处理**：采用优化的 PLINK 文件批量读取和索引技术
- ✅ **并行计算**：支持多线程并行执行，大幅提升分析速度

---

## 🆚 与 C++ 原版的对比

| 特性 | C++ SMR (v1.4.0) | pySMR (v0.3.0) | 说明 |
|------|------------------|----------------|------|
| **核心分析** | | | |
| 单 SNP SMR | ✅ | ✅ | 完全一致 |
| 多 SNP SMR | ✅ | ✅ *(新增)* | 算法完全对齐 |
| HEIDI 测试 | ✅ | ✅ | 完全一致 |
| **数据处理** | | | |
| 等位基因匹配 | 严格模式 | **鲁棒模式** | pySMR 支持链翻转/互补，数据利用率更高 |
| 输入格式 | BESD 二进制 | 灵活文本格式 | pySMR 更易用 |
| **高级功能** | | | |
| 并行计算 | ✅ (部分功能) | ✅ (全流程) | pySMR 将并行机制扩展到分析全流程 |
| **可视化** | | | |
| 绘图类型 | Locus / Effect | **+ Manhattan** | pySMR 提供更多图表 |
| 输出格式 | 仅 plot | **PDF / PNG** | 支持矢量图导出 |
| 基因标注 | ENSG ID | **Gene Symbol** | 支持自定义 ID 转换 |

---

## 🚀 快速开始

### 环境要求

- Python 3.8 或更高版本
- 操作系统：Linux、macOS、Windows

### 安装步骤

```bash
# 1. 克隆仓库
git clone https://github.com/CangmingDing/pySMR_Preview.git
cd pySMR_Preview

# 2. 安装依赖
pip install -r requirements.txt
```

### 基础用法

```bash
python main.py \
  --gwas "path/to/gwas_summary.txt" \
  --eqtl "path/to/eqtl_summary.csv" \
  --bfile "path/to/reference/1000G_EUR" \
  --out "output_prefix" \
  --threads 4
```

### 主要参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--gwas` | GWAS 摘要统计数据文件路径 | - |
| `--eqtl` | eQTL 摘要统计数据文件路径 | - |
| `--bfile` | PLINK 格式参考面板前缀（需 `.bed`, `.bim`, `.fam`） | - |
| `--ld-dir` | 分染色体 PLINK 文件目录（与 `--bfile` 二选一） | - |
| `--out` | 输出文件前缀 | - |
| `--threads` | CPU 线程数 | 1 |
| `--plot` | 是否生成可视化图表 | False |
| `--peqtl-smr` | 工具变量选择 P 值阈值（用于 Top-SNP 和 Multi-SNP） | 5e-8 |
| `--peqtl-heidi` | HEIDI 测试 SNP 选择 P 值阈值 | 1.57e-3 |
| `--cis-wind` | Cis 窗口大小（Kb） | 2000 |
| `--diff-freq` | 允许的最大等位基因频率差异 | 0.2 |
| `--no-strict-freq` | 禁用严格频率检查（不推荐） | False |
| **可视化参数** | | |
| `--plot-pdf` | 保存图表为 PDF 矢量格式 | False |
| `--add-gene-symbol` | 在图表中把 ID 转换为 Gene Symbol | False |
| `--gene-annotation` | ID 映射文件路径 (CSV/TSV, cols: ID, Symbol) | - |

---

## 📊 输出文件

分析完成后将生成以下文件：

- **`{out}.smr`**：包含所有探针 SMR 和 HEIDI 测试结果的制表符分隔文件
  - 包含字段：`probeID`, `topSNP`, `b_SMR`, `p_SMR`, **`p_SMR_multi`** *(新增)*, `p_HEIDI`, `nsnp_HEIDI`, 等
- **`{out}_plots/`**：包含显著基因的 Locus 图和 Effect Size 图的目录（需设置 `--plot`）

---

## 🎯 多 SNP SMR 测试详解

**什么是多 SNP SMR？**

传统的 SMR 测试仅使用区域内 P 值最小的单个 SNP（Top-SNP）作为工具变量。多 SNP SMR 则整合该基因所有显著独立的 eQTL 信号，通过以下步骤进行联合检验：

1. **选择候选 SNPs**：筛选 P 值小于阈值的所有 eQTL SNPs
2. **LD 修剪**：去除高度连锁的冗余 SNPs（r² > 0.95）
3. **计算统计量**：对每个独立 SNP 计算 SMR χ² 值并求和
4. **校正相关性**：使用 LD 矩阵特征值加权（Satterthwaite 近似）

**优势**：
- 📈 **更高的统计功效**：整合多个弱信号
- 🎯 **更稳健的结果**：减少单点估计的偶然性
- 🔬 **更精细的因果推断**：与 HEIDI 结合判断多效性

**输出字段**：
- `p_SMR_multi`：多 SNP 联合检验的 P 值
- `nsnp_multi`：用于多 SNP 测试的独立 SNP 数量

---

## 🔧 开发计划

- [x] ~~单 SNP SMR 测试~~
- [x] ~~HEIDI 异质性检验~~
- [x] ~~多 SNP SMR 测试~~ ✅ *已完成*
- [ ] 增强可视化效果（曼哈顿图、多层 Locus 图）
- [ ] 批处理模式优化
- [ ] Web 接口支持

---

## 📜 开源协议

本项目采用 [GPL-2.0 许可证](LICENSE)。

---

## 🙏 致谢

特别感谢杨剑教授及其团队开发的 [原始 SMR 软件](https://yanglab.westlake.edu.cn/software/smr/#Overview)，本项目严格遵循开源政策，在此基础上进行了 Python 实现并增强了数据处理鲁棒性。

---

## 📧 联系方式

- **作者**：Cangming
- **邮箱**：202201230726@163.com
- **微信**：CangMing-03

如有问题或建议，欢迎通过 [Issues](https://github.com/CangmingDing/pySMR_Preview/issues) 反馈！
