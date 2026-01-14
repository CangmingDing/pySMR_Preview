# pySMR_Preview Tutorial

This tutorial guides you through running an SMR analysis using `pySMR_Preview` (v0.3.0).

---

## Data Preparation

You need three key inputs:

1.  **GWAS Summary Statistics**: A text file with columns for SNP ID, Alleles (A1, A2), Beta, SE, P-value, and Frequency.
    *   *Example headers*: `SNP`, `A1`, `A2`, `b`, `se`, `p`, `freq`
2.  **eQTL Summary Statistics**: A text or CSV file for eQTLs, including Gene/Probe ID and position.
    *   *Example headers*: `SNP`, `Gene`, `A1`, `A2`, `beta`, `se`, `p_value`, `chrom`, `pos`
3.  **LD Reference Panel**: PLINK format (`.bed`, `.bim`, `.fam`). Usually 1000 Genomes data.

---

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

---

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

---

## Step 3: Understanding Multi-SNP SMR Results

**New in v0.3.0**: The output now includes Multi-SNP SMR test results!

### What is Multi-SNP SMR?

Traditional SMR uses only the **Top-SNP** (SNP with smallest eQTL P-value) as the instrument. Multi-SNP SMR integrates **multiple independent eQTL signals** in the region for a more powerful joint test.

**How it works**:
1. Selects all eQTL SNPs with `P < --peqtl-smr`
2. Performs LD pruning (removes SNPs with r² > 0.95 to ensure independence)
3. Computes a weighted chi-square sum statistic
4. Uses Satterthwaite approximation to calculate P-value

**Why use it?**
- ✅ **Higher Power**: Captures multiple weak signals
- ✅ **More Robust**: Less dependent on single-SNP noise
- ✅ **Better Inference**: Combined with HEIDI, provides stronger causal evidence

### Output Columns

The `.smr` output file now includes:

| Column | Description |
|--------|-------------|
| `p_SMR` | Single-SNP (Top-SNP) SMR P-value |
| **`p_SMR_multi`** | **Multi-SNP SMR P-value (NEW!)** |
| **`nsnp_multi`** | **Number of SNPs used in Multi-SNP test (NEW!)** |
| `p_HEIDI` | HEIDI heterogeneity test P-value |
| `nsnp_HEIDI` | Number of SNPs used in HEIDI test |

**Interpretation Guide**:
- If `p_SMR < 0.05` **AND** `p_SMR_multi < 0.05` **AND** `p_HEIDI > 0.01`:  
  → **Strong evidence** for causal effect (both single and multi-SNP tests agree)
- If `p_SMR_multi` is significant but `p_SMR` is not:  
  → Multi-SNP test captures weak distributed signals; consider region-level effect
- If `p_HEIDI < 0.01`:  
  → Heterogeneity detected; likely linkage or horizontal pleiotropy

---

## Understanding Plots

The tool generates two plots for each significantly associated gene:

### 1. Locus Plot (`*_locus.png`)

*   **Top Panel**: GWAS Manhattan plot for the region.
*   **Bottom Panel**: eQTL Manhattan plot.
*   **Diamond Marker**: The Top-SNP used for the Single-SNP SMR test.
*   **Gene Track**: Shows genes in the region with transcriptional direction arrows (→ or ←).

**What to look for**:
- Overlapping peaks in GWAS and eQTL plots suggest shared signal
- Multiple distinct peaks may indicate independent causal variants

### 2. Effect Size Plot (`*_effect.png`)

*   **Scatter plot** of GWAS Effect Size (y-axis) vs eQTL Effect Size (x-axis).
*   **Point Colors**: Represent LD ($r^2$) with the Top-SNP (red = high LD, blue = low LD).
*   **Regression Line**: The slope represents the SMR estimate ($b_{SMR}$).

**What to look for**:
- Points should cluster around the regression line (good SMR fit)
- High-LD SNPs (red) should align well; scattered points suggest heterogeneity
- Compare with `p_HEIDI`: if HEIDI fails and points are scattered, be cautious

### 3. SMR Manhattan Plot (`*_manhattan.png`)
*   **Genome-wide Overview**: Shows SMR P-values across all chromosomes.
*   **Shape Encoding**: Up-arrow ($\beta_{SMR} > 0$) vs Down-arrow ($\beta_{SMR} < 0$).
*   **Gene Labels**: Significant genes are automatically labeled.

---

## Output Quality & Customization

### PDF Vector Output
To generate publication-quality vector graphics (PDF) instead of PNG:

```bash
python main.py ... --plot --plot-pdf
```

### Gene Symbol Conversion
By default, plots use the IDs found in your eQTL file (often ENSG IDs). To display Gene Symbols:

1. Prepare a mapping file (CSV or TSV) with two columns: `[ID, Symbol]`.
2. Run with `--add-gene-symbol` and provide the file:

```bash
python main.py ... \
  --plot \
  --add-gene-symbol \
  --gene-annotation "human_gene_mapping.csv"
```


## Advanced Usage

### Custom Frequency Checking

By default, pySMR performs strict allele frequency checks (SNPs with freq difference > 0.2 are excluded). You can adjust this:

```bash
python main.py \
  --gwas "data.txt" \
  --eqtl "eqtl.csv" \
  --bfile "ref" \
  --out "out" \
  --diff-freq 0.3 \           # Allow larger freq difference
  --no-strict-freq            # Disable strict mode (NOT RECOMMENDED)
```

⚠️ **Warning**: Disabling `--no-strict-freq` may retain SNPs with strand errors or population mismatches.

### Chromosome-Specific Analysis

If your LD reference is split by chromosome:

```bash
python main.py \
  --gwas "data.txt" \
  --eqtl "eqtl.csv" \
  --ld-dir "./ld_ref/"  \      # Expects files like 1000G.EUR.1.bed, 1000G.EUR.2.bed, ...
  --out "out"
```

---

## Troubleshooting

**Q: Why are my Multi-SNP results `NA`?**  
A: This happens when:
- No SNPs pass the `--peqtl-smr` threshold
- Only 1 SNP remains after LD pruning (need ≥2 for multi-SNP test)
- LD matrix computation failed (no LD reference available)

**Q: My analysis is very slow!**  
A: Use `--threads` to enable parallel processing. For genome-wide data, 8-16 threads is optimal.

**Q: I get different results from C++ SMR**  
A: pySMR uses **more robust allele matching** (handles strand flips). This means:
- pySMR may retain more SNPs than C++ (if your data has strand mismatches)
- Core P-values should be very similar for well-formatted data
- If results diverge significantly, check your input data formatting

---

# pySMR_Preview 简易教程

本教程将指导您如何使用 `pySMR_Preview` (v0.3.0) 进行 SMR 分析。

---

## 数据准备

您需要准备三个关键输入文件：

1.  **GWAS 摘要统计数据**：包含 SNP ID、等位基因 (A1, A2)、Beta值、标准误 (SE)、P值和频率的文本文件。
    *   *常见表头示例*: `SNP`, `A1`, `A2`, `b`, `se`, `p`, `freq`
2.  **eQTL 摘要统计数据**：包含基因/探针 ID 和位置信息的 eQTL 数据。
    *   *常见表头示例*: `SNP`, `Gene`, `A1`, `A2`, `beta`, `se`, `p_value`, `chrom`, `pos`
3.  **LD 参考面板**：PLINK 格式文件 (`.bed`, `.bim`, `.fam`)。通常使用千人基因组 (1000 Genomes) 数据。

---

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

---

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

---

## 第三步：理解多 SNP SMR 结果

**v0.3.0 新增功能**：输出文件现在包含多 SNP SMR 测试结果！

### 什么是多 SNP SMR？

传统 SMR 仅使用**Top-SNP**（eQTL P 值最小的 SNP）作为工具变量。多 SNP SMR 则整合该区域内的**多个独立 eQTL 信号**进行更强大的联合检验。

**工作原理**：
1. 选择所有满足 `P < --peqtl-smr` 的 eQTL SNPs
2. 执行 LD 修剪（去除 r² > 0.95 的冗余 SNPs，确保独立性）
3. 计算加权卡方和统计量
4. 使用 Satterthwaite 近似计算 P 值

**为什么使用它？**
- ✅ **更高的统计功效**：捕获多个微弱信号
- ✅ **更稳健**：降低对单个 SNP 噪声的依赖
- ✅ **更好的推断**：结合 HEIDI，提供更强的因果证据

### 输出列说明

`.smr` 输出文件现在包含：

| 列名 | 说明 |
|------|------|
| `p_SMR` | 单 SNP (Top-SNP) SMR P 值 |
| **`p_SMR_multi`** | **多 SNP SMR P 值（新增！）** |
| **`nsnp_multi`** | **用于多 SNP 测试的 SNP 数量（新增！）** |
| `p_HEIDI` | HEIDI 异质性检验 P 值 |
| `nsnp_HEIDI` | 用于 HEIDI 测试的 SNP 数量 |

**解读指南**：
- 如果 `p_SMR < 0.05` **且** `p_SMR_multi < 0.05` **且** `p_HEIDI > 0.01`：  
  → **强因果证据**（单 SNP 和多 SNP 测试结果一致）
- 如果 `p_SMR_multi` 显著但 `p_SMR` 不显著：  
  → 多 SNP 测试捕获了分散的弱信号；考虑区域级效应
- 如果 `p_HEIDI < 0.01`：  
  → 检测到异质性；可能是连锁或水平多效性

---

## 图表解读

对于每个显著关联的基因，工具会生成两张图：

### 1. Locus 轨迹图 (`*_locus.png`)

*   **上图**: 该区域的 GWAS 曼哈顿图。
*   **下图**: eQTL 曼哈顿图。
*   **菱形标记**: 单 SNP SMR 测试使用的 Top-SNP。
*   **基因轨道**: 显示该区域内的基因及其转录方向箭头（→ 或 ←）。

**观察要点**：
- GWAS 和 eQTL 图中的峰值重叠表明共享信号
- 多个不同的峰值可能表明独立的因果变异

### 2. 效应值图 (`*_effect.png`)

*   **散点图**：GWAS 效应值（y轴）对 eQTL 效应值（x轴）。
*   **点的颜色**：表示与 Top-SNP 的 LD（$r^2$）强度（红色 = 高 LD，蓝色 = 低 LD）。
*   **回归线**：斜率代表 SMR 估计值（$b_{SMR}$）。

**观察要点**：
- 点应聚集在回归线周围（良好的 SMR 拟合）
- 高 LD SNPs（红色）应对齐良好；分散的点表明异质性
- 与 `p_HEIDI` 比较：如果 HEIDI 失败且点分散，需谨慎解释

### 3. SMR 曼哈顿图 (`*_manhattan.png`)
*   **全基因组概览**: 展示所有染色体上的 SMR P 值分布。
*   **形状编码**: 向上箭头 ($\beta_{SMR} > 0$) vs 向下箭头 ($\beta_{SMR} < 0$)。
*   **基因标注**: 显著基因会自动标注名称列表。

---

## 绘图自定义与发表级输出

### PDF 矢量输出
如果您需要用于文章发表的高清矢量图 (PDF)，而不是位图 (PNG)：

```bash
python main.py ... --plot --plot-pdf
```

### Gene Symbol 转换
默认情况下，图表使用 eQTL 文件中的 ID（通常是 ENSG ID）。如果您想显示基因符号 (Gene Symbol)：

1. 准备一个映射文件 (CSV 或 TSV)，包含两列：`[ID, Symbol]`。
2. 运行命令时添加 `--add-gene-symbol` 并指定文件路径：

```bash
python main.py ... \
  --plot \
  --add-gene-symbol \
  --gene-annotation "human_gene_mapping.csv"
```


## 高级用法

### 自定义频率检查

默认情况下，pySMR 执行严格的等位基因频率检查（频率差 > 0.2 的 SNPs 会被排除）。您可以调整：

```bash
python main.py \
  --gwas "data.txt" \
  --eqtl "eqtl.csv" \
  --bfile "ref" \
  --out "out" \
  --diff-freq 0.3 \           # 允许更大的频率差异
  --no-strict-freq            # 禁用严格模式（不推荐）
```

⚠️ **警告**：禁用 `--no-strict-freq` 可能会保留链错误或群体不匹配的 SNPs。

### 分染色体分析

如果您的 LD 参考面板按染色体分割：

```bash
python main.py \
  --gwas "data.txt" \
  --eqtl "eqtl.csv" \
  --ld-dir "./ld_ref/"  \      # 期望文件如 1000G.EUR.1.bed, 1000G.EUR.2.bed, ...
  --out "out"
```

---

## 常见问题

**问：为什么我的多 SNP 结果是 `NA`？**  
答：这发生在以下情况：
- 没有 SNP 通过 `--peqtl-smr` 阈值
- LD 修剪后只剩 1 个 SNP（多 SNP 测试需要 ≥2 个）
- LD 矩阵计算失败（无 LD 参考可用）

**问：我的分析很慢！**  
答：使用 `--threads` 启用并行处理。对于全基因组数据，8-16 个线程最优。

**问：我的结果与 C++ SMR 不同**  
答：pySMR 使用**更鲁棒的等位基因匹配**（处理链翻转）。这意味着：
- 如果您的数据有链不匹配，pySMR 可能保留比 C++ 更多的 SNPs
- 对于格式良好的数据，核心 P 值应该非常相似
- 如果结果显著不同，请检查输入数据格式
