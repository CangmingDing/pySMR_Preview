import pandas as pd
import numpy as np
import os
import argparse
from smrpy.plot import plot_locus, plot_effect_size
from smrpy.io import load_gwas, load_eqtl

def test_plot():
    # Load Real Data Snippets (for authenticity)
    print("Loading data for plotting test...")
    
    # We will just reuse the loading logic for a specific region to save time
    # Target Gene: ENSG00000132906.18 (Chr1)
    target_gene = "ENSG00000132906.18"
    target_chrom = 1
    
    gwas_path = "/Users/dingtunan/生信分析/肾癌/1-1网药/5基因组相关分析/0GWAS/肾癌1.txt"
    eqtl_path = "/Users/dingtunan/生信分析/肾癌/1-1网药/5基因组相关分析/0GWAS/kidney_eqtl.csv"
    
    # Load eQTL
    eqtl = pd.read_csv(eqtl_path) # Raw load to see structure
    # Use io.py logic
    from smrpy.io import standardize_columns
    eqtl = standardize_columns(eqtl)
    if 'tss_distance' in eqtl.columns and 'POS' in eqtl.columns:
        eqtl['TSS'] = eqtl['POS'] - eqtl['tss_distance']
        
    # Filter for target gene
    gene_data = eqtl[eqtl['GENE'] == target_gene].copy()
    if gene_data.empty:
        print(f"Gene {target_gene} not found in eQTL!")
        return

    print(f"Found {len(gene_data)} SNPs for {target_gene}")
    
    # Load GWAS for these SNPs
    gwas = pd.read_csv(gwas_path, sep='\t')
    gwas = standardize_columns(gwas)
    
    # Merge
    merged = pd.merge(gene_data, gwas, on='SNP', suffixes=('_EQTL', '_GWAS'))
    
    # Align (Simple for plotting test)
    # Assume A1/A2 match for now or just plot raw
    merged['B_EQTL_ALIGNED'] = merged['B_EQTL'] # Mock alignment
    
    top_snp_idx = merged['P_EQTL'].idxmin()
    top_row = merged.loc[top_snp_idx]
    
    # Calculate REAL SMR beta so the line passes through the red triangle
    # b_SMR = b_GWAS / b_eQTL
    real_b_smr = top_row['B_GWAS'] / top_row['B_EQTL_ALIGNED']
    
    mock_res = {
        'probeID': target_gene,
        'ProbeChr': target_chrom,
        'Probe_bp': gene_data['POS'].mean(),
        'topSNP': top_row['SNP'],
        'p_SMR': 3e-6,
        'b_SMR': real_b_smr, 
        'p_HEIDI': 0.05
    }
    
    # Region Genes (Mocking the track)
    region_genes = pd.DataFrame({
        'gene_id': [target_gene, 'NearbyGene1'],
        'Probe_bp': [merged['POS'].mean(), merged['POS'].mean() + 50000],
        'strand': ['+', '-']
    })
    
    out_dir = "debug_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    print("Generating Locus Plot...")
    plot_locus(merged, mock_res, f"{out_dir}/{target_gene}_locus.png", region_genes=region_genes)
    
    print("Generating Effect Plot...")
    # Mock LD vector for testing color gradient
    # Top SNP has r=1.0. Others have random r between 0 and 1.
    n_snps = len(merged)
    # Ensure Top SNP is 1.0
    mock_r_vec = np.random.uniform(0, 0.95, n_snps)
    # We need to find index of top snp in merged['SNP']
    top_snp_name = mock_res['topSNP']
    # But plot_effect_size expects r vector aligned to snp_list argument
    snp_list_arg = merged['SNP'].values
    
    # Set top SNP r to 1.0 in the mock vector (find index)
    top_idx = np.where(snp_list_arg == top_snp_name)[0]
    if len(top_idx) > 0:
        mock_r_vec[top_idx[0]] = 1.0
        
    plot_effect_size(merged, mock_res, f"{out_dir}/{target_gene}_effect.png", 
                     ld_r_vector=mock_r_vec, snp_list=snp_list_arg)
    
    print(f"Done! Check {out_dir}/")

if __name__ == "__main__":
    test_plot()
