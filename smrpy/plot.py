import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
from matplotlib.patches import Rectangle, Arrow
from matplotlib.collections import PatchCollection

import json

# --- Config Loader ---
def load_plot_config():
    """Load styling config from json, fallback to defaults if missing."""
    config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', 'plot_style.json')
    defaults = {
        "figure": {"figsize_locus": [11, 9], "figsize_effect": [7, 7], "figsize_manhattan": [14, 6], "figsize_qq": [6, 6], "figsize_volcano": [8, 6], "figsize_heidi": [8, 6], "figsize_multi": [8, 6], "dpi": 300},
        "fonts": {"family": "sans-serif", "sans_serif": ["Helvetica", "Arial"], "size_title": 14, "size_label": 11, "size_tick": 9},
        "colors": {
            "locus_gwas": "#777777", "locus_eqtl": "#A62A55", "locus_gene": "#8B8000", "locus_top_snp": "#4B0082", "locus_smr_line": "#D02090",
            "effect_errorbar": "#483D8B", "effect_regression": "#FF8C00", "effect_gradient": ["#4169E1", "#000080", "#101010"],
            "manhattan_palette": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"],
            "manhattan_sig_up": "#E41A1C", "manhattan_sig_down": "#377EB8", "manhattan_nonsig_alpha": 0.6,
            "qq_line": "red", "qq_ci": "grey", "volcano_sig": "#E41A1C", "volcano_nonsig": "grey", "heidi_pass": "#377EB8", "heidi_fail": "grey"
        },
        "markers": {"size_default": 20, "size_top": 80, "size_gene_arrow": 0.02}
    }
    
    if os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = json.load(f)
                # Simple merge (deep merge preferred but simple enough for now)
                for k, v in user_config.items():
                    if k in defaults and isinstance(v, dict):
                        defaults[k].update(v)
        except Exception as e:
            print(f"Warning: Failed to load plot config: {e}")
            
    return defaults

# Load Global Config once
PLOT_CONFIG = load_plot_config()

# --- Font Configuration ---
import matplotlib
matplotlib.rcParams['font.family'] = PLOT_CONFIG['fonts']['family']
if PLOT_CONFIG['fonts'].get('sans_serif'):
    matplotlib.rcParams['font.sans-serif'] = PLOT_CONFIG['fonts']['sans_serif'] + matplotlib.rcParams['font.sans-serif']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def convert_symbols(names, mapping=None):
    """
    Convert list of IDs to Symbols using mapping dict.
    If mapping is None, return as is (or strip version if needed).
    """
    if mapping is None:
        return names
    return [mapping.get(n, n) for n in names]

def plot_manhattan(smr_results, out_path, p_thresh=5e-8, mapping=None):
    """
    Generate SMR Manhattan Plot.
    Points colored by Chromosome.
    Shape indicates direction of effect (Beta > 0: Up, Beta < 0: Down).
    Significant genes are labeled.
    """
    if plt is None: return

    df = smr_results.copy()
    
    # Pre-process columns
    df['p_SMR'] = pd.to_numeric(df['p_SMR'], errors='coerce')
    df['b_SMR'] = pd.to_numeric(df['b_SMR'], errors='coerce')
    df['ProbeChr'] = pd.to_numeric(df['ProbeChr'], errors='coerce')
    df['Probe_bp'] = pd.to_numeric(df['Probe_bp'], errors='coerce')
    
    df = df.dropna(subset=['p_SMR', 'ProbeChr', 'Probe_bp'])
    df['logP'] = -np.log10(df['p_SMR'].replace(0, 1e-300))
    
    # Sort
    df = df.sort_values(['ProbeChr', 'Probe_bp'])
    
    # Cumulative Position Calculation
    chroms = sorted(df['ProbeChr'].unique())
    cum_pos_map = {}
    current_offset = 0
    ticks = []
    tick_labels = []
    
    df['rel_pos'] = df['Probe_bp']
    df['cum_pos'] = 0.0
    
    for chrom in chroms:
        idx = df['ProbeChr'] == chrom
        max_bp = df.loc[idx, 'Probe_bp'].max()
        min_bp = df.loc[idx, 'Probe_bp'].min()
        
        # Shift
        df.loc[idx, 'cum_pos'] = df.loc[idx, 'Probe_bp'] + current_offset
        
        # Tick in middle
        center_pos = current_offset + (max_bp + min_bp) / 2
        ticks.append(center_pos)
        tick_labels.append(str(int(chrom)))
        
        current_offset += max_bp
        
    # Plot setup
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_manhattan'])
    
    # Colors
    palette = PLOT_CONFIG['colors']['manhattan_palette']
    
    # Plot non-significant points first
    non_sig = df[df['p_SMR'] >= p_thresh]
    
    for i, chrom in enumerate(chroms):
        color = palette[i % len(palette)]
        subset = non_sig[non_sig['ProbeChr'] == chrom]
        if subset.empty: continue
        
        # Shape by beta direction
        # Matplotlib scatter is cleaner if we just plot dots for background
        # R script uses shape for ALL points. We can do that but it's slow for standard summary.
        # Let's use dots for non-sig to save rendering time/clutter, focus shapes on sig?
        # User prompt implies "参考 SMRr/R/plot.R". R script plots ALL with shapes.
        # Python scatter with varying markers is slow. We iterate by direction.
        
        # Up (Beta > 0)
        up = subset[subset['b_SMR'] > 0]
        if not up.empty:
            ax.scatter(up['cum_pos'], up['logP'], c=color, marker='^', s=15, alpha=0.6, edgecolors='none')
            
        # Down (Beta <= 0)
        down = subset[subset['b_SMR'] <= 0]
        if not down.empty:
            ax.scatter(down['cum_pos'], down['logP'], c=color, marker='v', s=15, alpha=0.6, edgecolors='none')

    # Plot Significant points
    sig = df[df['p_SMR'] < p_thresh]
    if not sig.empty:
        # Convert IDs here if mapping exists
        if mapping:
            sig['label'] = convert_symbols(sig['probeID'].values, mapping)
        else:
            sig['label'] = sig['probeID']

        ax.scatter(sig.loc[sig['b_SMR']>0, 'cum_pos'], sig.loc[sig['b_SMR']>0, 'logP'], 
                   c='red', marker='^', s=40, edgecolors='black', linewidth=0.5, zorder=5)
        ax.scatter(sig.loc[sig['b_SMR']<=0, 'cum_pos'], sig.loc[sig['b_SMR']<=0, 'logP'], 
                   c='red', marker='v', s=40, edgecolors='black', linewidth=0.5, zorder=5)
                   
        # Labels (using crude collision avoidance or logic)
        # Repel is hard in plain mpl without library. We'll do basic annotation for top hits.
        # Or limit to top 20 to avoid mess.
        
        top_n = sig.sort_values('p_SMR').head(30)
        texts = []
        for _, row in top_n.iterrows():
            t = ax.text(row['cum_pos'], row['logP'], row['label'], fontsize=8, fontstyle='italic', fontweight='bold')
            texts.append(t)
            
        # Optional: try adjustText if installed
        try:
            from adjustText import adjust_text
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        except ImportError:
            pass
            
    # Threshold Line
    ax.axhline(-np.log10(p_thresh), color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    # Styling
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=8, rotation=0)
    
    # Alternate shading background (optional, but R script doesn't seem to have it, just colors)
    
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P_{SMR})$")
    ax.set_title("SMR Manhattan Plot", fontweight='bold')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()


def plot_qq(smr_results, out_path):
    """
    Generate QQ Plot for SMR P-values.
    """
    if plt is None: return

    p_vals = pd.to_numeric(smr_results['p_SMR'], errors='coerce').dropna().values
    p_vals = p_vals[p_vals > 0] # Remove 0 or negative
    p_vals = np.sort(p_vals)
    n = len(p_vals)
    if n == 0: return

    # Theoretical Quantiles (Uniform)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(p_vals)
    
    # Sort Observed (smallest P -> largest -log10)
    observed = np.sort(observed)[::-1] # wait, smallest P is largest -log10.
    # standard: sort P small to large
    # expected: 1/(n+1) ... n/(n+1) -> small to large
    # observed should be sorted small to large p-values
    
    # Let's resort to be safe
    p_vals_sorted = np.sort(p_vals) # smallest first
    observed = -np.log10(p_vals_sorted)
    
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_qq'])
    
    ax.scatter(expected, observed, c='black', s=15, alpha=0.8)
    
    # Diagonal line
    max_val = max(np.max(expected), np.max(observed))
    ax.plot([0, max_val], [0, max_val], color=PLOT_CONFIG['colors']['qq_line'], linestyle='--')
    
    # CI (Optional, simplified 95%)
    # from scipy.stats import beta
    # For speed we skip complex CI shading unless requested.
    
    ax.set_xlabel(r'Expected $-\log_{10}(P)$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_ylabel(r'Observed $-\log_{10}(P)$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_title('SMR QQ Plot', fontweight=PLOT_CONFIG['fonts']['weight_title'])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()

def plot_volcano(smr_results, out_path, p_thresh=5e-8, mapping=None):
    """
    SMR Volcano Plot: b_SMR vs -log10(p_SMR)
    """
    if plt is None: return
    
    df = smr_results.copy()
    df['p_SMR'] = pd.to_numeric(df['p_SMR'], errors='coerce')
    df['b_SMR'] = pd.to_numeric(df['b_SMR'], errors='coerce')
    df = df.dropna(subset=['p_SMR', 'b_SMR'])
    
    df['logP'] = -np.log10(df['p_SMR'].replace(0, 1e-300))
    df['is_sig'] = df['p_SMR'] < p_thresh
    
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_volcano'])
    
    # Non-sig
    ax.scatter(df.loc[~df['is_sig'], 'b_SMR'], df.loc[~df['is_sig'], 'logP'],
               c=PLOT_CONFIG['colors']['volcano_nonsig'], s=10, alpha=0.5)
               
    # Sig
    sig_df = df[df['is_sig']]
    if not sig_df.empty:
        ax.scatter(sig_df['b_SMR'], sig_df['logP'],
                   c=PLOT_CONFIG['colors']['volcano_sig'], s=20, alpha=0.8)
                   
        # Label top hits
        top_n = sig_df.sort_values('p_SMR').head(10)
        
        if mapping: # Apply mapping for labels
             top_n['label'] = convert_symbols(top_n['probeID'].values, mapping)
        else:
             top_n['label'] = top_n['probeID']

        texts = []
        for _, row in top_n.iterrows():
            texts.append(ax.text(row['b_SMR'], row['logP'], row['label'], fontsize=8))
            
        try:
            from adjustText import adjust_text
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        except: pass

    ax.axhline(-np.log10(p_thresh), color=PLOT_CONFIG['colors']['volcano_sig'], linestyle='--', linewidth=0.8)
    ax.axvline(0, color='black', linestyle=':', linewidth=0.8)
    
    ax.set_xlabel(r'SMR Effect Size ($\beta_{SMR}$)', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_ylabel(r'$-\log_{10}(P_{SMR})$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_title('SMR Volcano Plot', fontweight=PLOT_CONFIG['fonts']['weight_title'])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()

def plot_smr_vs_heidi(smr_results, out_path, smr_thresh=5e-8, heidi_thresh=0.01):
    """
    SMR vs HEIDI Filtering Plot.
    X: -log10(P_SMR)
    Y: P_HEIDI (Linear scale, since we look for > 0.01)
       Or -log10(P_HEIDI) but usually we want high P_HEIDI.
    Let's plot P_HEIDI directly on Y (0 to 1). Or -log10 but inverted?
    Standard: -log10(P_SMR) vs P_HEIDI.
    Target: Bottom Right (High SMR Sig, High HEIDI P i.e. not significant heterogeneity)
    """
    if plt is None: return

    df = smr_results.copy()
    df['p_SMR'] = pd.to_numeric(df['p_SMR'], errors='coerce')
    df['p_HEIDI'] = pd.to_numeric(df['p_HEIDI'], errors='coerce')
    df = df.dropna(subset=['p_SMR', 'p_HEIDI'])
    
    df['logP_SMR'] = -np.log10(df['p_SMR'].replace(0, 1e-300))
    
    # Classify
    # Pass: SMR Sig AND HEIDI Not Sig ( > thresh)
    df['pass'] = (df['p_SMR'] < smr_thresh) & (df['p_HEIDI'] > heidi_thresh)
    # Fail HEIDI: SMR Sig AND HEIDI Sig ( < thresh)
    df['fail_heidi'] = (df['p_SMR'] < smr_thresh) & (df['p_HEIDI'] <= heidi_thresh)
    # Not Sig SMR
    df['nonsig_smr'] = df['p_SMR'] >= smr_thresh
    
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_heidi'])
    
    # Plot backgrounds
    ax.scatter(df.loc[df['nonsig_smr'], 'logP_SMR'], df.loc[df['nonsig_smr'], 'p_HEIDI'],
               c='lightgrey', s=10, alpha=0.3, label='Non-sig SMR')
               
    ax.scatter(df.loc[df['fail_heidi'], 'logP_SMR'], df.loc[df['fail_heidi'], 'p_HEIDI'],
               c=PLOT_CONFIG['colors']['heidi_fail'], marker='x', s=20, alpha=0.7, label='Heterogeneity (HEIDI Fail)')
               
    ax.scatter(df.loc[df['pass'], 'logP_SMR'], df.loc[df['pass'], 'p_HEIDI'],
               c=PLOT_CONFIG['colors']['heidi_pass'], marker='o', s=25, alpha=0.9, label='Robust Candidates (Pass)')
    
    # Threshold lines
    ax.axvline(-np.log10(smr_thresh), color='red', linestyle='--', linewidth=0.8, label='SMR Threshold')
    ax.axhline(heidi_thresh, color='blue', linestyle='--', linewidth=0.8, label='HEIDI Threshold')
    
    ax.set_xlabel(r'$-\log_{10}(P_{SMR})$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_ylabel(r'$P_{HEIDI}$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_title('SMR vs HEIDI Filtering', fontweight=PLOT_CONFIG['fonts']['weight_title'])
    
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()

def plot_single_vs_multi(smr_results, out_path):
    """
    Compare Single-SNP vs Multi-SNP P-values.
    Scatter plot: -log10(P_Single) vs -log10(P_Multi)
    """
    if plt is None: return

    df = smr_results.copy()
    if 'p_SMR_multi' not in df.columns: return
    
    df['p_SMR'] = pd.to_numeric(df['p_SMR'], errors='coerce')
    df['p_SMR_multi'] = pd.to_numeric(df['p_SMR_multi'], errors='coerce')
    df = df.dropna(subset=['p_SMR', 'p_SMR_multi'])
    
    df['logP_Single'] = -np.log10(df['p_SMR'].replace(0, 1e-300))
    df['logP_Multi'] = -np.log10(df['p_SMR_multi'].replace(0, 1e-300))
    
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_multi'])
    
    ax.scatter(df['logP_Single'], df['logP_Multi'], c='purple', alpha=0.5, s=15)
    
    # Diagonal y=x
    max_val = max(df['logP_Single'].max(), df['logP_Multi'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    ax.set_xlabel(r'Single-SNP $-\log_{10}(P)$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_ylabel(r'Multi-SNP $-\log_{10}(P)$', fontsize=PLOT_CONFIG['fonts']['size_label'])
    ax.set_title('Single-SNP vs Multi-SNP Power Comparison', fontweight=PLOT_CONFIG['fonts']['weight_title'])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()

def plot_locus(gene_data, smr_result, out_path, region_genes=None, mapping=None):
    """
    Generate an enhanced SMR Locus Plot (Triple Track).
    Top: GWAS -log10(P)
    Middle: eQTL -log10(P)
    Bottom: Gene Track with arrows
    """
    if plt is None:
        return

    # Determine Position Column
    pos_col = 'POS'
    if pos_col not in gene_data.columns:
        if 'POS_GWAS' in gene_data.columns: pos_col = 'POS_GWAS'
        elif 'POS_EQTL' in gene_data.columns: pos_col = 'POS_EQTL'
        else: return

    # Extract info
    top_snp = smr_result['topSNP']
    chrom = smr_result['ProbeChr']
    probe_id = smr_result['probeID']
    p_smr = smr_result['p_SMR']
    probe_bp = smr_result['Probe_bp']
    
    # Filter and prep data
    plot_data = gene_data.copy()
    plot_data['logP_GWAS'] = -np.log10(plot_data['P_GWAS'].replace(0, 1e-300))
    plot_data['logP_EQTL'] = -np.log10(plot_data['P_EQTL'].replace(0, 1e-300))
    
    # X-axis cleanup
    x_mb = plot_data[pos_col] / 1e6
    x_min, x_max = x_mb.min(), x_mb.max()
    pad = (x_max - x_min) * 0.05
    x_range = [x_min - pad, x_max + pad]

    # Style definitions from Config
    color_gwas = PLOT_CONFIG['colors']['locus_gwas']
    color_eqtl = PLOT_CONFIG['colors']['locus_eqtl']
    color_top = PLOT_CONFIG['colors']['locus_top_snp']
    color_smr_line = PLOT_CONFIG['colors']['locus_smr_line']
    color_gene = PLOT_CONFIG['colors']['locus_gene']
    
    # Fonts
    size_title = PLOT_CONFIG['fonts']['size_title']
    size_label = PLOT_CONFIG['fonts']['size_label']

    # Setup Figure
    fig = plt.figure(figsize=PLOT_CONFIG['figure']['figsize_locus'])
    gs = gridspec.GridSpec(3, 1, height_ratios=[4, 4, 2.5], hspace=0.1)
    
    # --- 1. GWAS Plot (Top) ---
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x_mb, plot_data['logP_GWAS'], 
                c=color_gwas, alpha=0.6, s=18, edgecolors='none', zorder=2)
    
    # Top SNP
    top_row = gene_data[gene_data['SNP'] == top_snp]
    if not top_row.empty:
        ax0.scatter(top_row[pos_col].iloc[0]/1e6, -np.log10(top_row['P_GWAS'].iloc[0]), 
                    c='none', marker='D', s=80, edgecolors=color_top, linewidths=1.5, zorder=5)
        ax0.axvline(x=top_row[pos_col].iloc[0]/1e6, color='blue', linestyle=':', alpha=0.3, linewidth=1)
    
    # SMR horizontal line
    ax0.axhline(y=-np.log10(p_smr) if p_smr > 0 else 0, color=color_smr_line, linestyle='--', alpha=0.7, linewidth=1)
    ax0.text(x_max, -np.log10(p_smr) + 0.3, f"$p_{{SMR}} = {p_smr:.1e}$", 
             ha='right', va='bottom', color=color_smr_line, fontsize=10, weight='bold')

    ax0.set_ylabel(r'$-\log_{10}(P\ \text{GWAS or SMR})$', fontsize=size_label)
    ax0.set_xlim(x_range)
    ax0.tick_params(labelbottom=False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    
    # --- 2. eQTL Plot (Middle) ---
    ax1 = plt.subplot(gs[1], sharex=ax0)
    ax1.scatter(x_mb, plot_data['logP_EQTL'], 
                color=color_eqtl, marker='x', alpha=0.7, s=20)
    
    if not top_row.empty:
         ax1.scatter(top_row[pos_col].iloc[0]/1e6, -np.log10(top_row['P_EQTL'].iloc[0]), 
                    c='none', marker='D', s=80, edgecolors=color_top, linewidths=1.5, zorder=5)

    ax1.set_ylabel(r'$-\log_{10}(P\ \text{eQTL})$', fontsize=size_label)
    ax1.tick_params(labelbottom=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # --- 3. Gene Track (Bottom) ---
    ax2 = plt.subplot(gs[2], sharex=ax0)
    
    if region_genes is not None and not region_genes.empty:
        unique_genes = region_genes.drop_duplicates(subset=['gene_id'])
        y_levels = [0.2, 0.6, 1.0, 1.4, 1.8]
        
        # Filter genes to be within Plot Range
        region_genes = region_genes[
            (region_genes['Probe_bp'] / 1e6 >= x_range[0]) & 
            (region_genes['Probe_bp'] / 1e6 <= x_range[1])
        ]

        for i, (idx, row) in enumerate(region_genes.iterrows()):
            g_name = row.get('gene_name', row['gene_id']).split('.')[0]
            g_pos = row.get('TSS', row.get('Probe_bp')) / 1e6
            
            # Strand logic: default to None/Unknown if not present or valid
            strand_raw = row.get('strand', '.')
            if strand_raw in ['+', '1', 1, 'fwd', 'forward']:
                direction = 1
            elif strand_raw in ['-', '-1', -1, 'rev', 'reverse']:
                direction = -1
            else:
                direction = 0 # No direction

            strand = row.get('strand', '+') # Keep for safe fallback if needed somewhere else? No, use logic above.
            
            y_lvl = y_levels[i % len(y_levels)]
            
            # Draw Arrow or Marker
            arrow_w = (x_max - x_min) * PLOT_CONFIG['markers']['size_gene_arrow']
            
            if direction != 0:
                ax2.arrow(g_pos, y_lvl, direction * arrow_w, 0, 
                          head_width=0.15, head_length=arrow_w*0.3, 
                          fc=color_gene, ec=color_gene, alpha=0.9, length_includes_head=True)
                lbl_x = g_pos + direction * arrow_w * 0.5
            else:
                # Draw a simple rectangle/marker for unknown strand
                rect_w = arrow_w * 0.8
                ax2.add_patch(Rectangle((g_pos - rect_w/2, y_lvl - 0.05), rect_w, 0.1, 
                                        fc=color_gene, ec=color_gene, alpha=0.6))
                lbl_x = g_pos 
            
            # Label
            ax2.text(lbl_x, y_lvl + 0.15, g_name, 
                     ha='center', va='bottom', fontsize=PLOT_CONFIG['fonts']['size_tick'], style='italic', color='#333333')
        
        ax2.set_ylim(0, 2.2)
    else:
        ax2.text(probe_bp/1e6, 0.5, f"Probe: {probe_id}", ha='center', va='center', fontsize=size_label, style='italic')
        ax2.set_ylim(0, 1)

    ax2.set_xlabel(f'Position on Chromosome {chrom} (Mb)', fontsize=size_label)
    ax2.get_yaxis().set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    # Main Title
    # Title
    display_id = probe_id
    if mapping:
         display_id = mapping.get(probe_id, probe_id)
         
    plt.suptitle(f"SMR Multi-omics Locus Plot: {display_id}", fontsize=size_title, y=0.98, weight=PLOT_CONFIG['fonts']['weight_title'])
    
    # Layout adjustment
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.08, right=0.95, hspace=0.15)
    # plt.tight_layout(rect=[0, 0.05, 1, 0.95]) # Caused UserWarning with GridSpec
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'])
    plt.close()

def plot_effect_size(gene_data, smr_result, out_path, ld_r_vector=None, snp_list=None, mapping=None):
    """
    Generate Effect Size Plot (Beta GWAS vs Beta eQTL).
    Matches SMR 'Effect scale' plot with error bars and r2 coloring.
    """
    if plt is None:
        return

    top_snp = smr_result['topSNP']
    b_smr = smr_result['b_SMR']
    
    # Extract aligned betas
    gene_data = gene_data.copy()
    
    fig, ax = plt.subplots(figsize=PLOT_CONFIG['figure']['figsize_effect'])
    
    # --- Style Config ---
    # Reference style: Clean, black ticks, no grid (or very subtle), ticks out
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=1, colors='black')
    
    # LD Coloring logic
    # Reference uses: Deep Blue (low r2) -> Black (high r2)
    # Let's create a custom colormap or use 'Blues' heavily truncated?
    # Actually the reference looks like: Blue -> Dark Blue -> Black
    import matplotlib.colors as mcolors
    # Custom cmap akin to the image
    cmap_colors = PLOT_CONFIG['colors']['effect_gradient'] # RoyalBlue -> Navy -> Near Black
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_blue_black", cmap_colors)

    snp_list_data = gene_data['SNP'].values
    r2_vals = np.zeros(len(gene_data)) # Default 0
    
    if ld_r_vector is not None and snp_list is not None:
        # Create mapping
        r2_map = dict(zip(snp_list, ld_r_vector**2))
        r2_vals = np.array([r2_map.get(s, 0.0) for s in snp_list_data])
        # Handle Nans
        r2_vals = np.nan_to_num(r2_vals)

    # Plot Non-Top SNPs
    # We want hollow circles with error bars crossing through
    
    # Separate Top SNP for special drawing
    is_top = gene_data['SNP'] == top_snp
    non_top = gene_data[~is_top]
    
    # Draw Error bars FIRST (zorder lower)
    # Color of error bars: nice slate blue/gray matching the points
    # We use a loop or ax.errorbar. 
    # To get the color mapping for error bars matching points is hard with single call if we want gradient.
    # Actually reference has UNIFORM error bar color (blueish-grey) for non-top points?
    # Looking closely at reference: Error bars seem same color as the point (Blue).
    # Easier: Uniform styling for error bars, say 'navy' or '#333333'
    
    # Scatter plot for colorbar mapping
    sc = ax.scatter(non_top['B_EQTL_ALIGNED'], non_top['B_GWAS'], 
                    c=r2_vals[~is_top], cmap=custom_cmap, vmin=0, vmax=1,
                    s=50, alpha=0.9, marker='o', edgecolors='face', linewidths=1.5, zorder=5)
    
    # We need "hollow" circles? The reference is solid? 
    # PROMPT SAYS: "普通点：空心圆圈（带叉号/误差线）" -> User perception.
    # Looking at the artifact: The points look like circles with a cross inside.
    # They look like Open Circles with a cross. 
    # Matplotlib marker for that is NOT standard.
    # Or just `facecolors='none', edgecolors=color`.
    # Let's try facecolors='none' but we need to map edgecolors to value.
    
    # Re-draw scatter manually to handle 'hollow' with mapped edge color
    norm = mcolors.Normalize(0, 1)
    mapped_colors = custom_cmap(norm(r2_vals[~is_top]))
    
    ax.scatter(non_top['B_EQTL_ALIGNED'], non_top['B_GWAS'], 
               c='none', s=60, marker='o', linewidths=1.5, edgecolors=mapped_colors, zorder=5)
               
    # Add cross inside? Or just error bars effectively form the cross?
    # Typically regular error bars + open circle looks like the reference.
    
    # Error Bars for non-top
    # We color error bars same as the edge color? Or uniform?
    # Reference: Error bars are same blue-ish color.
    # Let's iterate efficiently? No, too slow. Use uniform 'navy' equivalent color.
    # Or simply:
    ax.errorbar(non_top['B_EQTL_ALIGNED'], non_top['B_GWAS'], 
                xerr=non_top['SE_EQTL'], yerr=non_top['SE_GWAS'], 
                fmt='none', ecolor=PLOT_CONFIG['colors']['effect_errorbar'], alpha=0.3, elinewidth=0.8, capsize=0, zorder=4)

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, aspect=20)
    cbar.set_label(r'$r^2$', rotation=0, labelpad=15, y=0.5)
    cbar.set_ticks([0.14, 0.31, 0.48, 0.66, 0.83, 1.00]) # Mimic reference or auto
    # Let's use auto but cleaner
    
    # Highlight Top SNP
    # Red hollow triangle with cross
    top_row = gene_data[is_top]
    if not top_row.empty:
        # Error bars
        ax.errorbar(top_row['B_EQTL_ALIGNED'], top_row['B_GWAS'], 
                    xerr=top_row['SE_EQTL'], yerr=top_row['SE_GWAS'], 
                    fmt='none', ecolor='red', alpha=0.8, elinewidth=1.2, capsize=0, zorder=9)
        # Marker: Hollow red triangle
        ax.scatter(top_row['B_EQTL_ALIGNED'], top_row['B_GWAS'], 
                   c='none', s=120, marker='^', linewidths=2.0, edgecolors='red', zorder=10, label='top cis-eQTL')

    # Regression Line (SMR slope)
    # Orange dashed
    xlim = ax.get_xlim()
    # Expand slightly to fit line nicely
    extra = (xlim[1] - xlim[0]) * 0.1
    line_x = np.linspace(xlim[0] - extra, xlim[1] + extra, 100)
    line_y = b_smr * line_x
    ax.plot(line_x, line_y, '--', color=PLOT_CONFIG['colors']['effect_regression'], linewidth=2, alpha=0.8, zorder=3, label='SMR Estimate') # Deep orange
    
    # dummy legend entries to match style
    # Top SNP
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='^', color='w', label='top cis-eQTL',
               markerfacecolor='none', markeredgecolor='red', markersize=10, markeredgewidth=1.5),
        Line2D([0], [0], marker='o', color='w', label='cis-eQTL',
               markerfacecolor='none', markeredgecolor='#000080', markersize=8, markeredgewidth=1.5)
    ]
    ax.legend(handles=legend_elements, loc='best', frameon=False)

    ax.set_xlabel('eQTL effect sizes', fontsize=PLOT_CONFIG['fonts']['size_label'], labelpad=10)
    ax.set_ylabel('GWAS effect sizes', fontsize=PLOT_CONFIG['fonts']['size_label'], labelpad=10)
    
    # Title: Probe (Gene)
    title_str = f"{smr_result['probeID']}" 
    if mapping:
        title_str = mapping.get(smr_result['probeID'], title_str)
    
    # If gene name available separately, append it. Usually probeID is gene in our pipe.
    ax.set_title(title_str, fontsize=PLOT_CONFIG['fonts']['size_label'], loc='left', pad=15)
    
    ax.grid(False) # Reference has no grid

    plt.tight_layout()
    plt.savefig(out_path, dpi=PLOT_CONFIG['figure']['dpi'], bbox_inches='tight')
    plt.close()

