# pySMR_Preview

<div align="center">

**基于 Python 的SMR分析工具（预览版）**

[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)

[English Version](README_EN.md) | 中文文档

</div>

---

## 📖 项目简介

**pySMR_Preview** 是对杨剑教授团队开发的原始 [SMR 软件](https://yanglab.westlake.edu.cn/software/smr/#Overview) 的 Python 实现（预览版）。本项目严格遵循开源政策，目前已复刻 SMR 和 HEIDI 测试的核心功能。

> **⚠️ 重要说明**  
> 当前版本为**预览版**，仅实现了基础的 SMR 和 HEIDI 分析流程。完整的功能复刻（包括高级可视化、批处理优化等）将在后续版本中持续更新。

本工具旨在为整合 GWAS 和 eQTL 摘要数据提供灵活、高效的 Python 解决方案。

---

## ✨ 核心功能

- ✅ **SMR 测试**：使用单个遗传变异作为工具变量，检测基因表达与复杂性状之间的关联
- ✅ **HEIDI 测试**：检测依赖工具变量中的异质性，以区分多效性（Pleiotropy）和连锁（Linkage）
- ✅ **基础可视化**：生成 Locus 轨迹图和效应值（Effect Size）散点图
- ✅ **高效 LD 处理**：采用优化的 PLINK 文件批量读取和索引技术
- ✅ **并行计算**：支持多线程并行执行，大幅提升分析速度

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
| `--peqtl-smr` | 工具变量选择 P 值阈值 | 5e-8 |
| `--cis-wind` | Cis 窗口大小（Kb） | 2000 |

---

## 📊 输出文件

分析完成后将生成以下文件：

- **`{out}.smr`**：包含所有探针 SMR 和 HEIDI 测试结果的制表符分隔文件
- **`{out}_plots/`**：包含显著基因的 Locus 图和 Effect Size 图的目录（需设置 `--plot`）

---

## 🔧 开发计划

- [ ] 增强可视化效果（曼哈顿图、多层 Locus 图）
- [ ] 批处理模式优化
- [ ] 更多质控参数
- [ ] Web 接口支持

---

## 📜 开源协议

本项目采用 [GPL-2.0 许可证](LICENSE)。

---

## 🙏 致谢

特别感谢杨剑教授及其团队开发的 [原始 SMR 软件](https://yanglab.westlake.edu.cn/software/smr/#Overview)，本项目严格遵循开源政策，在此基础上进行了 Python 实现。

---

## 📧 联系方式

- **作者**：Cangming
- **邮箱**：202201230726@bucm.edu.cn
- **微信**：CangMing-03

如有问题或建议，欢迎通过 [Issues](https://github.com/CangmingDing/pySMR_Preview/issues) 反馈！
