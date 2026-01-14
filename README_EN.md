# pySMR_Preview

<div align="center">

**Python Implementation of Summary-data-based Mendelian Randomization (v0.3.0 - Full Feature Release)**

[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)

[中文文档](README.md) | English Version

</div>

---

## 📖 About

**pySMR_Preview** is a Python implementation (v0.3.0) of the original [SMR software](https://yanglab.westlake.edu.cn/software/smr/#Overview) developed by Prof. Jian Yang's group. This project strictly adheres to open-source policies and has **fully replicated the core functionality** of SMR and HEIDI tests, with enhanced allele matching robustness.

> **✅ Feature Parity Status**  
> The current version achieves **complete alignment** with the C++ original at the analysis level:
> - ✅ Single-SNP SMR test (Top-SNP SMR)
> - ✅ **Multi-SNP SMR test (Set-based SMR)** - **New Feature!**
> - ✅ HEIDI heterogeneity test
> - ✅ Enhanced allele matching (supports strand flips and complementary base matching)

This tool aims to provide a flexible and efficient Python solution for integrating GWAS and eQTL summary data.

---

## ✨ Key Features

- ✅ **Single-SNP SMR Test**: Tests for association between gene expression and complex traits using a single genetic variant as an instrument
- ✅ **🎯 Multi-SNP SMR Test** *(New)*:  
  Joint testing using multiple independent eQTL signals (Set-based test), computing P-values via Satterthwaite approximation on weighted chi-square sums, significantly improving statistical power
- ✅ **HEIDI Test**: Tests for heterogeneity in dependent instruments to distinguish pleiotropy from linkage
- ✅ **Advanced Allele Matching**: Automatically handles strand flips and base complements (A/G ↔ T/C), maximizing data utilization
- ✅ **Complete Visualization**: Generates Locus plots, Effect Size scatter plots, and **Manhattan Plots**
- ✅ **Publication-Ready Plots**: Supports vector PDF output, Helvetica fonts, and **Gene Symbol conversion**
- ✅ **Efficient LD Handling**: Optimized batch reading and indexing for PLINK reference files
- ✅ **Parallel Processing**: Multi-threaded execution support for significant speedup

---

## 🆚 Comparison with C++ Original

| Feature | C++ SMR (v1.4.0) | pySMR (v0.3.0) | Notes |
|---------|------------------|----------------|-------|
| **Core Analysis** | | | |
| Single-SNP SMR | ✅ | ✅ | Fully aligned |
| Multi-SNP SMR | ✅ | ✅ *(New)* | Algorithm fully aligned |
| HEIDI Test | ✅ | ✅ | Fully aligned |
| **Data Processing** | | | |
| Allele Matching | Strict mode | **Robust mode** | pySMR supports strand flips/complements, higher data retention |
| Input Format | BESD binary | Flexible text | pySMR more user-friendly |
| **Advanced Features** | | | |
| Parallel Computing | ✅ (Partial) | ✅ (Full Workflow) | pySMR extends parallelism to global analysis |
| **Visualization** | | | |
| Plot Types | Locus / Effect | **+ Manhattan** | pySMR provides more plot types |
| Output Format | Plot only | **PDF / PNG** | Supports vector graphics export |
| Gene Labeling | ENSG ID | **Gene Symbol** | Supports custom ID conversion |

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
| `--peqtl-smr` | P-value threshold for instrument selection (Top-SNP & Multi-SNP) | 5e-8 |
| `--peqtl-heidi` | P-value threshold for HEIDI test SNP selection | 1.57e-3 |
| `--cis-wind` | Cis-window size (Kb) | 2000 |
| `--diff-freq` | Maximum allowed allele frequency difference | 0.2 |
| `--no-strict-freq` | Disable strict frequency check (NOT RECOMMENDED) | False |
| **Visualization Arguments** | | |
| `--plot-pdf` | Save plots as PDF vector graphics | False |
| `--add-gene-symbol` | Convert IDs to Gene Symbols in plots | False |
| `--gene-annotation` | Path to ID mapping file (CSV/TSV, cols: ID, Symbol) | - |

---

## 📊 Output Files

After analysis completion, the following files will be generated:

- **`{out}.smr`**: Tab-delimited file containing SMR and HEIDI test results for all probes
  - Fields include: `probeID`, `topSNP`, `b_SMR`, `p_SMR`, **`p_SMR_multi`** *(New)*, `p_HEIDI`, `nsnp_HEIDI`, etc.
- **`{out}_plots/`**: Directory containing Locus and Effect Size plots for significant genes (requires `--plot` flag)

---

## 🎯 Multi-SNP SMR Test Explained

**What is Multi-SNP SMR?**

Traditional SMR tests use only the single SNP with the smallest P-value (Top-SNP) in the region as an instrument. Multi-SNP SMR integrates all significant independent eQTL signals for joint testing through the following steps:

1. **Select Candidates**: Filter all eQTL SNPs with P-value below threshold
2. **LD Pruning**: Remove highly correlated redundant SNPs (r² > 0.95)
3. **Compute Statistic**: Calculate SMR χ² values for each independent SNP and sum them
4. **Correct for Correlation**: Weight by eigenvalues of LD matrix (Satterthwaite approximation)

**Advantages**:
- 📈 **Higher Statistical Power**: Integrates multiple weak signals
- 🎯 **More Robust Results**: Reduces randomness of single-point estimates
- 🔬 **Finer Causal Inference**: Combined with HEIDI to assess pleiotropy

**Output Fields**:
- `p_SMR_multi`: P-value of multi-SNP joint test
- `nsnp_multi`: Number of independent SNPs used in multi-SNP test

---

## 🔧 Development Roadmap

- [x] ~~Single-SNP SMR test~~
- [x] ~~HEIDI heterogeneity test~~
- [x] ~~Multi-SNP SMR test~~ ✅ *Completed*
- [ ] Enhanced visualizations (Manhattan plots, multi-track locus plots)
- [ ] Batch processing optimization
- [ ] Web interface support

---

## 📜 License

This project is licensed under the [GPL-2.0 License](LICENSE).

---

## 🙏 Acknowledgements

Special thanks to Prof. Jian Yang and his team for developing the [original SMR software](https://yanglab.westlake.edu.cn/software/smr/#Overview). This project strictly follows open-source policies, implements the methodology in Python, and enhances data processing robustness.

---

## 📧 Contact

- **Author**: Cangming
- **Email**: 202201230726@163.com
- **WeChat**: CangMing-03

For questions or suggestions, please submit via [Issues](https://github.com/CangmingDing/pySMR_Preview/issues)!
