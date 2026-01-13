import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
from matplotlib.patches import Rectangle, Arrow
from matplotlib.collections import PatchCollection

def plot_locus(gene_data, smr_result, out_path, region_genes=None):
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

    # Style definitions (matching SMR reference)
    color_gwas = '#777777' # Gray
    color_eqtl = '#A62A55' # Maroon/Pink for eQTL dots
    color_top = '#4B0082'  # Indigo/Purple for top SNP
    color_smr_line = '#D02090' # Pinkish line
    color_gene = '#8B8000' # Golden for arrows

    # Setup Figure
    fig = plt.figure(figsize=(11, 9))
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

    ax0.set_ylabel(r'$-\log_{10}(P\ \text{GWAS or SMR})$', fontsize=11)
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

    ax1.set_ylabel(r'$-\log_{10}(P\ \text{eQTL})$', fontsize=11)
    ax1.tick_params(labelbottom=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # --- 3. Gene Track (Bottom) ---
    ax2 = plt.subplot(gs[2], sharex=ax0)
    
    if region_genes is not None and not region_genes.empty:
        unique_genes = region_genes.drop_duplicates(subset=['gene_id'])
        y_levels = [0.2, 0.6, 1.0, 1.4, 1.8]
        
        for i, (idx, row) in enumerate(unique_genes.iterrows()):
            g_name = row.get('gene_name', row['gene_id']).split('.')[0]
            g_pos = row.get('TSS', row.get('Probe_bp')) / 1e6
            strand = row.get('strand', '+')
            
            y_lvl = y_levels[i % len(y_levels)]
            
            # Draw Arrow
            arrow_w = (x_max - x_min) * 0.015
            direction = 1 if strand in ['+', '1', 1] else -1
            ax2.arrow(g_pos, y_lvl, direction * arrow_w, 0, 
                      head_width=0.15, head_length=arrow_w*0.3, fc=color_gene, ec=color_gene, alpha=0.9)
            
            # Label
            ax2.text(g_pos + direction * arrow_w * 0.5, y_lvl + 0.1, g_name, 
                     ha='center', va='bottom', fontsize=9, style='italic', color='#333333')
        
        ax2.set_ylim(0, 2.2)
    else:
        ax2.text(probe_bp/1e6, 0.5, f"Probe: {probe_id}", ha='center', va='center', fontsize=10, style='italic')
        ax2.set_ylim(0, 1)

    ax2.set_xlabel(f'Position on Chromosome {chrom} (Mb)', fontsize=11)
    ax2.get_yaxis().set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    # Main Title
    plt.suptitle(f"SMR Multi-omics Locus Plot: {probe_id}", fontsize=14, y=0.98, weight='bold')
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(out_path, dpi=250)
    plt.close()

def plot_effect_size(gene_data, smr_result, out_path, ld_r_vector=None, snp_list=None):
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
    
    fig, ax = plt.subplots(figsize=(7, 7))
    
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
    cmap_colors = ["#4169E1", "#000080", "#101010"] # RoyalBlue -> Navy -> Near Black
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
                fmt='none', ecolor='#483D8B', alpha=0.3, elinewidth=0.8, capsize=0, zorder=4)

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
    ax.plot(line_x, line_y, '--', color='#FF8C00', linewidth=2, alpha=0.8, zorder=3, label='SMR Estimate') # Deep orange
    
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

    ax.set_xlabel('eQTL effect sizes', fontsize=12, labelpad=10)
    ax.set_ylabel('GWAS effect sizes', fontsize=12, labelpad=10)
    
    # Title: Probe (Gene)
    title_str = f"{smr_result['probeID']}" 
    # If gene name available separately, append it. Usually probeID is gene in our pipe.
    ax.set_title(title_str, fontsize=12, loc='left', pad=15)
    
    ax.grid(False) # Reference has no grid

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

