import argparse
import pandas as pd
import numpy as np
import os
from smrpy.io import load_gwas, load_eqtl, get_ld_matrix, get_ld_vector, read_plink1_bin, load_file_chunked, clear_snp_cache, _SNP_INDEX_CACHE
from smrpy.alignment import align_stats
from smrpy.stats import smr_test, heidi_test, smr_multi_snp_test
from smrpy.plot import plot_locus, plot_effect_size, plot_manhattan, plot_qq, plot_volcano, plot_smr_vs_heidi, plot_single_vs_multi
from tqdm import tqdm
import gc
from joblib import Parallel, delayed

# --- Helper Functions for Processing ---

# Global workers for parallel execution
_worker_G_ref = None
_worker_snp_map = None

def init_worker(G_ref, snp_map):
    """Initialize worker process with shared read-only data (if using process backend, likely not used for G due to size, relying on threading backend)"""
    global _worker_G_ref
    global _worker_snp_map
    _worker_G_ref = G_ref
    _worker_snp_map = snp_map

def process_gene(gene, chrom, gene_data, args, G_ref_chrom=None):
    """
    Process a single gene.
    If G_ref_chrom is not None, it is used (sequential).
    If None, tries to use global _worker_G_ref (parallel).
    """
    # Resolve reference
    if G_ref_chrom is None:
        global _worker_G_ref
        G_ref_chrom = _worker_G_ref
    
    # Ensure cache is populated if using global G
    # (Note: In threading backend, cache is shared. In process backend, it's local to worker)
    if G_ref_chrom is not None:
         obj_id = id(G_ref_chrom)
         if obj_id not in _SNP_INDEX_CACHE:
             # Try to restore map if we have it in global
             global _worker_snp_map
             if _worker_snp_map:
                 _SNP_INDEX_CACHE[obj_id] = (G_ref_chrom, _worker_snp_map)

    # Select Instruments (P < peqtl_smr)
    candidates = gene_data[gene_data['P_EQTL'] <= args.peqtl_smr]
    if candidates.empty:
        return None
        
    top_snp_idx = candidates['P_EQTL'].idxmin()
    top_snp_row = gene_data.loc[top_snp_idx]
    
    # SMR Test
    # Calculate aligned beta first
    b_xy, se_xy, p_smr = smr_test(top_snp_row['B_EQTL_ALIGNED'], 
                                  top_snp_row['SE_EQTL'],
                                  top_snp_row['B_GWAS'],
                                  top_snp_row['SE_GWAS'])
    
    # --- Multi-SNP SMR Test ---
    p_smr_multi = np.nan
    nsnp_multi = 0
    
    # 1. Select Candidates (P < peqtl_smr)
    # candidates df already defined above
    if args.smr_multi and len(candidates) >= 1 and G_ref_chrom is not None:
        try:
             # Get LD for all candidates
            snp_list_multi = candidates['SNP'].values
            ld_mat_multi, snps_in_ld_multi = get_ld_matrix(G_ref_chrom, snp_list_multi, chrom=chrom)
            
            if ld_mat_multi is not None and len(snps_in_ld_multi) > 0:
                # Align data to LD matrix order
                multi_data = candidates.set_index('SNP').loc[snps_in_ld_multi]
                
                # Pruning for stability (Remove exact duplicates/very high LD)
                # Greedy selection based on P-value
                # If we don't prune, close-to-singular matrices cause issues for eigen-decomposition
                # Using args.ld_multi_snp default 0.1
                
                prune_thresh = args.ld_multi_snp
                sorted_idx = np.argsort(multi_data['P_EQTL'].values)
                selected_indices = []
                
                # Keep track of 'removed' indices
                removed_mask = np.zeros(len(sorted_idx), dtype=bool)
                
                for i in sorted_idx:
                    if removed_mask[i]:
                        continue
                    
                    selected_indices.append(i)
                    
                    # Mark neighbors as removed
                    r2_vec = ld_mat_multi[i, :]**2
                    to_remove = np.where(r2_vec > prune_thresh)[0]
                    removed_mask[to_remove] = True
                    
                selected_indices = np.array(selected_indices)
                
                if len(selected_indices) > 0:
                     # Subset Data
                    b_zx_m = multi_data.iloc[selected_indices]['B_EQTL_ALIGNED'].values
                    se_zx_m = multi_data.iloc[selected_indices]['SE_EQTL'].values
                    b_zy_m = multi_data.iloc[selected_indices]['B_GWAS'].values
                    se_zy_m = multi_data.iloc[selected_indices]['SE_GWAS'].values
                    
                    ld_mat_final_m = ld_mat_multi[np.ix_(selected_indices, selected_indices)]
                    
                    p_smr_multi, nsnp_multi = smr_multi_snp_test(b_zx_m, se_zx_m, b_zy_m, se_zy_m, ld_mat_final_m)

        except Exception as e:
            # print(f"Multi-SMR Error: {e}")
            pass

    
    # HEIDI Test
    heidi_snps_data = gene_data[gene_data['P_EQTL'] < args.peqtl_heidi]
    p_heidi = np.nan
    nsnp_heidi = 0
    
    # LD Logic
    if len(heidi_snps_data) >= 3 and G_ref_chrom is not None:
        try:
            snp_list = heidi_snps_data['SNP'].values
            ld_mat, snps_in_ld = get_ld_matrix(G_ref_chrom, snp_list, chrom=chrom)
            
            if ld_mat is not None and len(snps_in_ld) >= 3:
                if top_snp_row['SNP'] in snps_in_ld:
                    top_idx_in_ld = np.where(snps_in_ld == top_snp_row['SNP'])[0][0]
                    r2_with_top = ld_mat[top_idx_in_ld]**2
                    
                    # Filter
                    mask = (r2_with_top < args.ld_upper_limit) & (r2_with_top > args.ld_lower_limit)
                    mask[top_idx_in_ld] = True 
                    
                    selected = np.where(mask)[0]
                    
                    # Max M Pruning
                    if len(selected) > args.heidi_max_m:
                        snps_final_candidates = snps_in_ld[selected]
                        sub_df = heidi_snps_data[heidi_snps_data['SNP'].isin(snps_final_candidates)]
                        sub_df = sub_df[sub_df['SNP'] != top_snp_row['SNP']]
                        sub_df = sub_df.sort_values('P_EQTL')
                        
                        keepers = np.append(sub_df.head(args.heidi_max_m - 1)['SNP'].values, top_snp_row['SNP'])
                        final_mask = np.isin(snps_in_ld, keepers)
                        selected = np.where(final_mask)[0]

                    if len(selected) >= 3:
                        ld_mat_final = ld_mat[np.ix_(selected, selected)]
                        snps_final = snps_in_ld[selected]
                        h_data = heidi_snps_data[heidi_snps_data['SNP'].isin(snps_final)].set_index('SNP').loc[snps_final]
                        
                        p_heidi, nsnp_heidi = heidi_test(h_data['B_EQTL_ALIGNED'].values,
                                                         h_data['SE_EQTL'].values,
                                                         h_data['B_GWAS'].values,
                                                         h_data['SE_GWAS'].values,
                                                         ld_mat_final)
        except Exception:
            pass
    
    res_entry = {
        'probeID': gene, 
        'ProbeChr': chrom,
        'Probe_bp': top_snp_row['POS'],
        'topSNP': top_snp_row['SNP'],
        'A1': top_snp_row['A1_GWAS'],
        'A2': top_snp_row['A2_GWAS'],
        'Freq': top_snp_row.get('FREQ_GWAS', np.nan),
        'b_GWAS': top_snp_row['B_GWAS'],
        'se_GWAS': top_snp_row['SE_GWAS'],
        'p_GWAS': top_snp_row['P_GWAS'],
        'b_eQTL': top_snp_row['B_EQTL_ALIGNED'],
        'se_eQTL': top_snp_row['SE_EQTL'],
        'p_eQTL': top_snp_row['P_EQTL'],
        'b_SMR': b_xy,
        'se_SMR': se_xy,
        'p_SMR': p_smr,
        'p_SMR_multi': p_smr_multi,
        'nsnp_multi': nsnp_multi,
        'p_HEIDI': p_heidi,
        'nsnp_HEIDI': nsnp_heidi
    }
    return res_entry

# --- Main Pipeline ---

def main():
    parser = argparse.ArgumentParser(description="Python implementation of SMR and HEIDI (v0.3.0).")
    
    # Input/Output
    parser.add_argument("--gwas", required=True, help="GWAS summary statistics")
    parser.add_argument("--eqtl", required=True, help="eQTL summary statistics")
    parser.add_argument("--ld-dir", help="Directory containing chromosome-separated PLINK files (1000G.EUR.{CHR})")
    parser.add_argument("--bfile", help="Single PLINK binary prefix (e.g. /path/to/ref if single file)")
    parser.add_argument("--out", required=True, help="Output file prefix")
    
    # Parallelization
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads to use (default 1 = sequential)")

    # Flags
    parser.add_argument("--low-memory", action='store_true', help="Use chunked reading to save RAM (slower)")
    
    # Thresholds
    parser.add_argument("--peqtl-smr", type=float, default=5e-8, help="P-value threshold for SMR (default 5e-8)")
    parser.add_argument("--peqtl-heidi", type=float, default=1.57e-3, help="P-value threshold for HEIDI SNPs (default 1.57e-3)")
    parser.add_argument("--maf", type=float, default=0.01, help="MAF filtering threshold (default 0.01)")
    parser.add_argument("--diff-freq", type=float, default=0.2, help="Maximum allele frequency difference allowed (default 0.2)")
    parser.add_argument("--cis-wind", type=int, default=2000, help="Cis-window in Kb (default 2000)")
    
    # Frequency Strictness (New)
    parser.add_argument("--no-strict-freq", action='store_true', help="Disable strict frequency check (NOT RECOMMENDED). If set, SNPs with freq difference > threshold will be kept if alleles match.")
    parser.add_argument("--freq-prop-tol", type=float, default=0.05, help="Tolerance of proportion of SNPs failing frequency check (default 0.05)")

    # HEIDI Params
    parser.add_argument("--heidi-min-m", type=int, default=3, help="Min SNPs (m) for HEIDI test (default 3)")
    parser.add_argument("--heidi-max-m", type=int, default=20, help="Max SNPs for HEIDI test")
    parser.add_argument("--ld-upper-limit", type=float, default=0.9, help="LD max r2 for HEIDI SNPs (default 0.9)")
    parser.add_argument("--ld-lower-limit", type=float, default=0.05, help="LD min r2 for HEIDI SNPs (default 0.05)")
    
    # Multi-SNP SMR Params (Parity with C++)
    parser.add_argument("--smr-multi", action='store_true', help="Enable Multi-SNP SMR test (default: False)")
    parser.add_argument("--ld-multi-snp", type=float, default=0.1, help="LD r2 threshold for Multi-SNP pruning (default 0.1 matches C++ SMR)")

    
    # Visualization
    parser.add_argument("--plot", action='store_true', help="Generate locus plots for significant results")
    parser.add_argument("--plot-threshold", type=float, default=5e-8, help="P-SMR threshold for plotting")
    parser.add_argument("--plot-pdf", action='store_true', help="Save plots as PDF instead of PNG")
    parser.add_argument("--add-gene-symbol", action='store_true', help="Convert probe IDs to Gene Symbols in plots")
    parser.add_argument("--gene-annotation", help="Path to CSV/TSV file for ID mapping (columns: probeID, geneSymbol)")

    args = parser.parse_args()

    # Derived arguments logic
    args.strict_freq = not args.no_strict_freq  # Default is True unless flag is set

    print("\n--- SMRpy v0.3.0 Pipeline Started ---")
    if args.threads > 1:
        print(f"Running in PARALLEL mode with {args.threads} threads.")
    else:
        print("Running in SEQUENTIAL mode.")
    
    if args.bfile and args.ld_dir:
        print("Warning: Both --bfile and --ld-dir provided. Using --bfile.")
    if not args.bfile and not args.ld_dir:
        print("Error: Must provide either --bfile or --ld-dir for LD reference.")
        return

    # --- Step 1: Load eQTL (Usually smaller) ---
    print(f"Loading eQTL: {args.eqtl}")
    eqtl = load_eqtl(args.eqtl)
    if 'tss_distance' in eqtl.columns and 'POS' in eqtl.columns:
        print(f"  Calculating TSS and identifying gene regions...")
        eqtl['TSS'] = eqtl['POS'] - eqtl['tss_distance']
    
    if 'tss_distance' in eqtl.columns:
        print(f"  Filtering eQTLs by cis-window ({args.cis_wind}Kb)...")
        limit = args.cis_wind * 1000
        n_before = len(eqtl)
        eqtl = eqtl[eqtl['tss_distance'].abs() <= limit]
        print(f"  Kept {len(eqtl)} / {n_before} probes.")
        
    # Get List of SNPs present in eQTL
    eqtl_snps = eqtl['SNP'].unique()
    print(f"  Unique eQTL SNPs: {len(eqtl_snps)}")
    
    # --- Step 2: Load GWAS (Usually huge) ---
    print(f"Loading GWAS: {args.gwas}")
    
    if args.low_memory:
        print("  [Low Memory Mode] filtering GWAS on-the-fly...")
        gwas = load_file_chunked(args.gwas, filter_col='SNP', filter_values=eqtl_snps)
    else:
        # Verbose False to keep logs clean
        gwas = load_gwas(args.gwas, verbose=False)
        gwas = gwas[gwas['SNP'].isin(eqtl_snps)]
    
    # Basic GWAS QC
    if 'FREQ' in gwas.columns:
        n_before = len(gwas)
        gwas = gwas[gwas['FREQ'].between(args.maf, 1 - args.maf)]
        print(f"  GWAS MAF Filter ({args.maf}): {n_before} -> {len(gwas)} SNPs after merge.")
        
    if gwas.empty:
        print("Error: No GWAS SNPs remained after filtering/matching.")
        return
    
    # --- Step 3: Alignment ---
    print("Aligning GWAS and eQTL...")
    data = align_stats(gwas, eqtl, freq_thresh=args.diff_freq, strict_freq=args.strict_freq)
    
    if data.empty:
        print("No overlapping SNPs found after alignment and QC.")
        return
        
    print(f"Total valid SNP-Gene pairs Analysis set: {len(data)}")

    # --- Step 4: Analysis ---
    
    single_ref = None
    snp_map_chrom = None # Global map cache
    
    if args.bfile:
        print(f"Loading Global LD Reference: {args.bfile}")
        try:
            single_ref = read_plink1_bin(args.bfile + ".bed", 
                                   args.bfile + ".bim", 
                                   args.bfile + ".fam", verbose=False)
            
            from smrpy.io import load_bim_with_progress, _SNP_INDEX_CACHE
            snp_map_chrom = load_bim_with_progress(args.bfile + ".bim")
            obj_id = id(single_ref)
            _SNP_INDEX_CACHE[obj_id] = (single_ref, snp_map_chrom)
            
        except Exception as e:
            print(f"Error loading PLINK: {e}")
            return

    chroms = sorted(data['CHR'].unique())
    all_results = []
    
    # Plotting setup
    plot_ext = "pdf" if args.plot_pdf else "png"
    gene_mapping = None
    
    if args.add_gene_symbol:
        if args.gene_annotation:
            try:
                print(f"Loading gene annotation: {args.gene_annotation}")
                # Try simple load: assumes header
                anno_df = pd.read_csv(args.gene_annotation, sep=None, engine='python')
                # Assume first two columns are ID and Symbol if not named specifically
                if 'probeID' in anno_df.columns and 'geneSymbol' in anno_df.columns:
                    gene_mapping = dict(zip(anno_df['probeID'], anno_df['geneSymbol']))
                else:
                    # Fallback to index 0 and 1
                    gene_mapping = dict(zip(anno_df.iloc[:, 0], anno_df.iloc[:, 1]))
                print(f"  Loaded {len(gene_mapping)} mappings.")
            except Exception as e:
                print(f"  Error loading annotation: {e}")
    
    if args.plot:
        plot_dir = f"{args.out}_plots"
        os.makedirs(plot_dir, exist_ok=True)

    for chrom in chroms:
        try:
            chrom_int = int(chrom)
        except:
            chrom_int = chrom
            
        print(f"\nProcessing Chromosome {chrom}...")
        
        # Subset data
        chr_data = data[data['CHR'] == chrom]
        if chr_data.empty:
            continue
            
        genes = chr_data['GENE'].unique()
        print(f"  - {len(genes)} genes.")
        
        G_ref_chrom = single_ref
        
        # Load LD if not global
        if G_ref_chrom is None:
             # Folder mode
            bfile_prefix = os.path.join(args.ld_dir, f"1000G.EUR.{chrom}")
            if not os.path.exists(bfile_prefix + ".bed"):
                print(f"  - LD reference not found ({bfile_prefix}), skipping...")
                continue
            try:
                # OPTIMIZED LOADING
                from smrpy.io import load_bim_with_progress, _SNP_INDEX_CACHE
                
                snp_map_chrom = load_bim_with_progress(bfile_prefix + ".bim")
                G_ref_chrom = read_plink1_bin(bfile_prefix + ".bed", 
                                       bfile_prefix + ".bim", 
                                       bfile_prefix + ".fam", verbose=False)
                
                obj_id = id(G_ref_chrom)
                _SNP_INDEX_CACHE[obj_id] = (G_ref_chrom, snp_map_chrom)
                
            except Exception as e:
                print(f"  - Error loading PLINK: {e}")
                continue

        # Execution Selection
        results_chrom = []
        
        if args.threads > 1:
            # Parallel Execution
            # Using sharedmem (threading) is crucial for G_ref sharing
            # We initialize workers with G_ref (although threading backend shares memory anyway, init logic is good practice)
            init_worker(G_ref_chrom, snp_map_chrom)
            
            gene_tasks = []
            for gene in genes:
                gene_data = chr_data[chr_data['GENE'] == gene]
                gene_tasks.append(delayed(process_gene)(gene, chrom, gene_data, args, None)) # Pass None for G to use global
                
            try:
                # Using 'threading' backend to share G_ref_chrom without pickling
                results_chrom = Parallel(n_jobs=args.threads, require='sharedmem')(tqdm(gene_tasks, desc=f"Chr {chrom} (Parallel)", unit="gene"))
                # Filter Nones
                results_chrom = [r for r in results_chrom if r is not None]
            except Exception as e:
                print(f"Parallel execution failed: {e}")
                continue
        else:
            # Sequential Execution
            for gene in tqdm(genes, desc=f"Chr {chrom} (Sequential)"):
                gene_data = chr_data[chr_data['GENE'] == gene]
                res = process_gene(gene, chrom, gene_data, args, G_ref_chrom=G_ref_chrom) # Pass G explicitly
                if res:
                    results_chrom.append(res)
        
        # Post-Processing: Aggregation & Plotting (Sequential Safe)
        for res_entry in results_chrom:
            all_results.append(res_entry)
            
            if args.plot and res_entry['p_SMR'] < args.plot_threshold:
                 # Re-retrieve plotting data (inefficient? yes but safe)
                gene = res_entry['probeID']
                gene_data = chr_data[chr_data['GENE'] == gene]
                
                # Get nearby genes
                gene_pos = res_entry['Probe_bp']
                window_bp = args.cis_wind * 1000
                region_eqtl = eqtl[(eqtl['CHR'] == chrom) & (eqtl['POS'].between(gene_pos - window_bp, gene_pos + window_bp))]
                
                if not region_eqtl.empty:
                    agg_dict = {
                        'TSS': 'first',
                        'POS': 'mean',
                    }
                    
                    # Try to capture strand if it exists
                    strand_col = None
                    for col in ['strand', 'STRAND', 'orientation', 'ori']:
                        if col in region_eqtl.columns:
                            strand_col = col
                            break
                    
                    if strand_col:
                         agg_dict[strand_col] = 'first'
                    
                    region_genes = region_eqtl.groupby('GENE').agg(agg_dict).reset_index().rename(columns={'GENE': 'gene_id'})
                    
                    if strand_col:
                        region_genes['strand'] = region_genes[strand_col]
                    else:
                        region_genes['strand'] = '.' # Unknown
                        
                    region_genes['Probe_bp'] = region_genes['POS']
                else:
                    region_genes = None

                locus_plot_name = f"{gene}_{chrom}_locus.{plot_ext}"
                locus_plot_path = os.path.join(plot_dir, locus_plot_name)
                plot_locus(gene_data, res_entry, locus_plot_path, region_genes=region_genes, mapping=gene_mapping)
                
                effect_plot_name = f"{gene}_{chrom}_effect.{plot_ext}"
                effect_plot_path = os.path.join(plot_dir, effect_plot_name)
                
                try:
                    all_snps = gene_data['SNP'].values
                    ld_vec, snp_list_ld = get_ld_vector(G_ref_chrom, res_entry['topSNP'], all_snps, chrom=chrom)
                    plot_effect_size(gene_data, res_entry, effect_plot_path, ld_r_vector=ld_vec, snp_list=snp_list_ld, mapping=gene_mapping)
                except Exception:
                    plot_effect_size(gene_data, res_entry, effect_plot_path, mapping=gene_mapping)

        # Cleanup
        clear_snp_cache()
        gc.collect()

    if all_results:
        out_df = pd.DataFrame(all_results)
        output_path = f"{args.out}.smr"
        out_df.to_csv(output_path, sep='\t', index=False)
        print(f"\nAnalysis complete. Results saved to {output_path}")
        if args.plot:
            # Generate Manhattan Plot
            print("Generating Manhattan Plot...")
            manhattan_path = os.path.join(plot_dir, f"{os.path.basename(args.out)}_manhattan.{plot_ext}")
            plot_manhattan(out_df, manhattan_path, p_thresh=args.plot_threshold, mapping=gene_mapping)
            
            # Additional Plots
            print("Generating additional diagnostic plots...")
            plot_qq(out_df, os.path.join(plot_dir, f"{os.path.basename(args.out)}_qq.{plot_ext}"))
            plot_volcano(out_df, os.path.join(plot_dir, f"{os.path.basename(args.out)}_volcano.{plot_ext}"), p_thresh=args.plot_threshold, mapping=gene_mapping)
            plot_smr_vs_heidi(out_df, os.path.join(plot_dir, f"{os.path.basename(args.out)}_smr_vs_heidi.{plot_ext}"), smr_thresh=args.plot_threshold)
            plot_single_vs_multi(out_df, os.path.join(plot_dir, f"{os.path.basename(args.out)}_single_vs_multi.{plot_ext}"))
            
            print(f"Plots saved to {plot_dir}/")
    else:
        print("\nNo significant results to save.")

if __name__ == "__main__":
    main()
