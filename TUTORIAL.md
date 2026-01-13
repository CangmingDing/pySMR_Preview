# pySMR_Preview Tutorial

This tutorial guides you through running an SMR analysis using `pySMR_Preview`.

## Data Preparation

You need three key inputs:

1.  **GWAS Summary Statistics**: A text file with columns for SNP ID, Alleles (A1, A2), Beta, SE, P-value, and Frequency.
    *   *Example headers*: `SNP`, `A1`, `A2`, `b`, `se`, `p`, `freq`
2.  **eQTL Summary Statistics**: A text or CSV file for eQTLs, including Gene/Probe ID and position.
    *   *Example headers*: `SNP`, `Gene`, `A1`, `A2`, `beta`, `se`, `p_value`, `chrom`, `pos`
3.  **LD Reference Panel**: PLINK format (`.bed`, `.bim`, `.fam`). Usually 1000 Genomes data.

## Step 1: Sequential Analysis (Test Run)

It is recommended to first run a sequential analysis (default) to ensure data loads correctly and parameters are reasonable.

```bash
python main.py \
  --gwas "./data/GWAS_Kidney.txt" \
  --eqtl "./data/eQTL_Kidney.csv" \
  --bfile "./reference/1000G_EUR" \
  --out "result_test" \
  --plot \
  --peqtl-smr 1e-5 \
  --cis-wind 500
```

*   **`--peqtl-smr 1e-5`**: We use a relaxed threshold for testing to ensure we find some candidate instruments.
*   **`--cis-wind 500`**: Defines the region around the gene to look for SNPs (500Kb).

## Step 2: High-Performance Parallel Analysis

Once the test run is successful, switch to parallel mode for the full analysis. This is critical for genome-wide data.

```bash
python main.py \
  --gwas "./data/GWAS_Kidney.txt" \
  --eqtl "./data/eQTL_Kidney.csv" \
  --bfile "./reference/1000G_EUR" \
  --out "result_final" \
  --plot \
  --peqtl-smr 5e-8 \
  --threads 8
```

*   **`--threads 8`**: This enables parallel processing on 8 cores. The script uses shared memory to efficiently access the large LD reference without crashing RAM.

## Understanding Plots

The tool generates two plots for each significant result:

1.  **Locus Plot (`*_locus.png`)**:
    *   **Top Panel**: GWAS Manhattan plot for the region.
    *   **Bottom Panel**: eQTL Manhattan plot.
    *   **Diamond**: The top SNP used for the SMR test.
    *   **Gene Track**: Shows genes in the region with transcriptional direction arrows.

2.  **Effect Size Plot (`*_effect.png`)**:
    *   Scatter plot of GWAS Effect Size (y-axis) vs eQTL Effect Size (x-axis).
    *   Points are colored by LD ($r^2$) with the top SNP.
    *   The slope represents the SMR estimate ($b_{SMR}$).

---

# pySMR_Preview 简易教程

本教程将指导您如何使用 `pySMR_Preview` 进行 SMR 分析。

## 数据准备

您需要准备三个关键输入文件：

1.  **GWAS 摘要统计数据**：包含 SNP ID、等位基因 (A1, A2)、Beta值、标准误 (SE)、P值和频率的文本文件。
    *   *常见表头示例*: `SNP`, `A1`, `A2`, `b`, `se`, `p`, `freq`
2.  **eQTL 摘要统计数据**：包含基因/探针 ID 和位置信息的 eQTL 数据。
    *   *常见表头示例*: `SNP`, `Gene`, `A1`, `A2`, `beta`, `se`, `p_value`, `chrom`, `pos`
3.  **LD 参考面板**：PLINK 格式文件 (`.bed`, `.bim`, `.fam`)。通常使用千人基因组 (1000 Genomes) 数据。

## 第一步：单线程测试运行

建议首先运行默认的单线程模式，以检查数据读取是否正确，参数设置是否合理。

```bash
python main.py \
  --gwas "./data/GWAS_Kidney.txt" \
  --eqtl "./data/eQTL_Kidney.csv" \
  --bfile "./reference/1000G_EUR" \
  --out "result_test" \
  --plot \
  --peqtl-smr 1e-5 \
  --cis-wind 500
```

*   **`--peqtl-smr 1e-5`**: 在测试时使用较宽松的阈值，确保能筛选到一些候选工具变量。
*   **`--cis-wind 500`**: 定义基因周围寻找 SNP 的范围 (500Kb)。

## 第二步：高性能并行分析

测试成功后，开启并行模式进行全基因组分析。这对于大规模数据至关重要。

```bash
python main.py \
  --gwas "./data/GWAS_Kidney.txt" \
  --eqtl "./data/eQTL_Kidney.csv" \
  --bfile "./reference/1000G_EUR" \
  --out "result_final" \
  --plot \
  --peqtl-smr 5e-8 \
  --threads 8
```

*   **`--threads 8`**: 启用 8 个线程并行处理。脚本利用共享内存技术，在多核读取大型 LD 矩阵时不会导致内存溢出。

## 图表解读

对于每个显著结果，工具会生成两张图：

1.  **Locus 轨迹图 (`*_locus.png`)**:
    *   **上图**: 该区域的 GWAS 曼哈顿图。
    *   **下图**: eQTL 曼哈顿图。
    *   **菱形点**: SMR 测试使用的 Top SNP。
    *   **基因轨道**: 显示该区域内的基因及其转录方向箭头。

2.  **效应值图 (`*_effect.png`)**:
    *   GWAS 效应值 (y轴) 对 eQTL 效应值 (x轴) 的散点图。
    *   点的颜色根据其与 Top SNP 的 LD ($r^2$) 强度显示。
    *   斜率代表 SMR 估计值 ($b_{SMR}$)。
