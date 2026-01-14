import pandas as pd
import numpy as np

def align_stats(gwas_df, eqtl_df, freq_thresh=0.2, strict_freq=True):
    """
    Align GWAS and eQTL summary statistics.
    Includes:
    1. Allele Flipping (consistent effect alleles)
    2. Strand Check (complementary alleles)
    3. Frequency Check (for palindromic SNPs and general QC)
    """
    # Merge on SNP ID
    merged = pd.merge(gwas_df, eqtl_df, on='SNP', suffixes=('_GWAS', '_EQTL'))
    
    if merged.empty:
        return merged

    # Helper for complement
    compl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    def get_compl(a):
        return compl.get(a, a)

    def check_alignment(row):
        # GWAS info
        a1g, a2g = row['A1_GWAS'], row['A2_GWAS']
        # eQTL info
        a1e, a2e = row['A1_EQTL'], row['A2_EQTL']
        
        # Frequency info (if available)
        freq_g = row.get('FREQ_GWAS', np.nan)
        freq_e = row.get('FREQ_EQTL', np.nan)
        
        # Check alignment
        # return: normalization_factor (1, -1, or 0 if remove)
        
        factor = 0
        
        # Case 1: Direct Match
        if a1g == a1e and a2g == a2e:
            factor = 1
        # Case 2: Flipped Match
        elif a1g == a2e and a2g == a1e:
            factor = -1
        # Case 3: Complement Match (Strand flip)
        elif a1g == get_compl(a1e) and a2g == get_compl(a2e):
            factor = 1
        # Case 4: Complement Flipped
        elif a1g == get_compl(a2e) and a2g == get_compl(a1e):
            factor = -1
            
        # Frequency Check mechanism
        # Only if frequencies are available in both datasets
        if factor != 0 and pd.notnull(freq_g) and pd.notnull(freq_e):
            
            freq_diff = 0
            if factor == 1:
                freq_diff = abs(freq_g - freq_e)
            else:
                freq_diff = abs(freq_g - (1.0 - freq_e))
                
            if freq_diff > freq_thresh:
                # SMR logic: Remove if frequency difference is too large. 
                # This indicates likely strand error or different populations.
                if strict_freq:
                    return 0 # Exclude
                # If strict_freq is False, we keep it (risk of strand error)
                
            # Ambiguous (Palindromic) Check logic from SMR (A/T or C/G)
            # If alleles are palindromic, we strictly rely on frequency.
            # If freq is close to 0.5, we can't tell strands apart -> Remove.
            # Usually SMR warns if maf > 0.4 for palindromic. 
            # We will implement the Diff Check primarily.
            
        return factor

    merged['ALIGN_FACTOR'] = merged.apply(check_alignment, axis=1)
    
    # Filter
    merged = merged[merged['ALIGN_FACTOR'] != 0].copy()
    
    # Adjust Beta
    merged['B_EQTL_ALIGNED'] = merged['B_EQTL'] * merged['ALIGN_FACTOR']
    
    return merged
