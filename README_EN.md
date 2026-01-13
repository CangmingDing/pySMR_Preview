# pySMR_Preview

<div align="center">

**Python Implementation of Summary-data-based Mendelian Randomization (Preview Version)**

[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)

[中文文档](README.md) | English Version

</div>

---

## 📖 About

**pySMR_Preview** is a Python implementation (preview version) of the original [SMR software](https://yanglab.westlake.edu.cn/software/smr/#Overview) developed by Prof. Jian Yang's group. This project strictly adheres to open-source policies and currently replicates the core functionality of SMR and HEIDI tests.

> **⚠️ Important Notice**  
> This is a **preview version** that implements only the fundamental SMR and HEIDI analysis pipeline. Full feature replication (including advanced visualizations, batch processing optimizations, etc.) will be continuously updated in future releases.

This tool aims to provide a flexible and efficient Python solution for integrating GWAS and eQTL summary data.

---

## ✨ Key Features

- ✅ **SMR Test**: Tests for association between gene expression and complex traits using a single genetic variant as an instrument
- ✅ **HEIDI Test**: Tests for heterogeneity in dependent instruments to distinguish pleiotropy from linkage
- ✅ **Basic Visualization**: Generates locus plots and effect size scatter plots
- ✅ **Efficient LD Handling**: Optimized batch reading and indexing for PLINK reference files
- ✅ **Parallel Processing**: Multi-threaded execution support for significant speedup

---

## 🚀 Quick Start

### Requirements

- Python 3.8 or higher
- OS: Linux, macOS, Windows

### Installation

```bash
# 1. Clone repository
git clone https://github.com/CangmingDing/pySMR_Preview.git
cd pySMR_Preview

# 2. Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```bash
python main.py \
  --gwas "path/to/gwas_summary.txt" \
  --eqtl "path/to/eqtl_summary.csv" \
  --bfile "path/to/reference/1000G_EUR" \
  --out "output_prefix" \
  --threads 4
```

### Key Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--gwas` | Path to GWAS summary statistics | - |
| `--eqtl` | Path to eQTL summary statistics | - |
| `--bfile` | PLINK binary prefix (requires `.bed`, `.bim`, `.fam`) | - |
| `--ld-dir` | Directory with chromosome-separated PLINK files (alternative to `--bfile`) | - |
| `--out` | Output file prefix | - |
| `--threads` | Number of CPU threads | 1 |
| `--plot` | Generate visualization plots | False |
| `--peqtl-smr` | P-value threshold for instrument selection | 5e-8 |
| `--cis-wind` | Cis-window size (Kb) | 2000 |

---

## 📊 Output Files

After analysis completion, the following files will be generated:

- **`{out}.smr`**: Tab-delimited file containing SMR and HEIDI test results for all probes
- **`{out}_plots/`**: Directory containing Locus and Effect Size plots for significant genes (requires `--plot` flag)

---

## 🔧 Development Roadmap

- [ ] Enhanced visualizations (Manhattan plots, multi-track locus plots)
- [ ] Batch processing optimization
- [ ] Additional QC parameters
- [ ] Web interface support

---

## 📜 License

This project is licensed under the [GPL-2.0 License](LICENSE).

---

## 🙏 Acknowledgements

Special thanks to Prof. Jian Yang and his team for developing the [original SMR software](https://yanglab.westlake.edu.cn/software/smr/#Overview). This project strictly follows open-source policies and implements the methodology in Python.

---

## 📧 Contact

- **Author**: Cangming
- **Email**: 202201230726@163.com
- **WeChat**: CangMing-03

For questions or suggestions, please submit via [Issues](https://github.com/CangmingDing/pySMR_Preview/issues)!
