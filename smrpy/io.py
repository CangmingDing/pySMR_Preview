import pandas as pd
import numpy as np
import os
from pandas_plink import read_plink1_bin
import re


def standardize_columns(df):
    """
    Auto-detect and rename columns to internal standard:
    SNP, A1, A2, B, SE, P, FREQ, CHR, POS, GENE
    """
    # Define mapping of standard names to possible user variations (regex)
    # ORDER MATTERS in the regex list for priority if multiple columns match same pattern? 
    # Actually we iterate columns.
    
    # We explicitly handle 'rsid' vs 'snp' priority logic separately below.
    
    maps = {
        'A1': [r'^a1$', r'^effect_allele$', r'^ref$', r'^allele1$'],
        'A2': [r'^a2$', r'^other_allele$', r'^alt$', r'^allele2$'],
        'B': [r'^b$', r'^beta$', r'^or$', r'^effect$', r'^slope$'],
        'SE': [r'^se$', r'^stderr$', r'^slope_se$'],
        'P': [r'^p$', r'^pval$', r'^p_val$', r'^pvalue$', r'^pval_nominal$'],
        'FREQ': [r'^freq$', r'^frq$', r'^maf$', r'^eaf$', r'^af$'],
        'CHR': [r'^chr$', r'^chromosome$', r'^chrom$'],
        'POS': [r'^pos$', r'^bp$', r'^position$', r'^base_pair$'],
        'GENE': [r'^gene$', r'^gene_id$', r'^probe$', r'^probe_id$']
    }
    
    cols = df.columns.tolist()
    rename_dict = {}
    
    # Special Handling for SNP ID priority
    # If 'rsid' exists, use it as SNP. 
    # If not, try 'snp', 'marker', etc.
    
    # Check for RSID first
    rsid_col = None
    for c in cols:
        if re.match(r'^rsid$', c, re.IGNORECASE):
            rsid_col = c
            break
            
    if rsid_col:
        rename_dict[rsid_col] = 'SNP'
    else:
        # Fallback to generic 'snp'
        for c in cols:
            if re.match(r'^snp$', c, re.IGNORECASE) or re.match(r'^marker$', c, re.IGNORECASE):
                rename_dict[c] = 'SNP'
                break

    # Map other columns
    for standard, patterns in maps.items():
        if standard in cols: continue 
        if standard == 'SNP' and 'SNP' in rename_dict.values(): continue
        
        found = False
        for pat in patterns:
            for c in cols:
                # Avoid renaming a column that is already mapped (e.g. rsid mapped to SNP)
                if c in rename_dict: continue
                
                if re.match(pat, c, re.IGNORECASE):
                    rename_dict[c] = standard
                    found = True
                    break
            if found: break
            
    if rename_dict:
        print(f"  Auto-mapped columns: {rename_dict}")
        df = df.rename(columns=rename_dict)
        
    return df

def load_gwas(path, verbose=True):
    if verbose:
        print(f"  Reading {path}...")
    try:
        with open(path, 'r') as f:
            line = f.readline()
            if '\t' in line: sep = '\t'
            elif ',' in line: sep = ','
            else: sep = r'\s+'
        df = pd.read_csv(path, sep=sep)
    except:
        df = pd.read_csv(path, sep=r'\s+')
        
    df = standardize_columns(df)
    return df

def load_eqtl(path, verbose=True):
    if verbose:
        print(f"  Reading {path}...")
    try:
        with open(path, 'r') as f:
            line = f.readline()
            if '\t' in line: sep = '\t'
            elif ',' in line: sep = ','
            else: sep = r'\s+'
        df = pd.read_csv(path, sep=sep)
    except:
        df = pd.read_csv(path, sep=r'\s+')
    

    df = standardize_columns(df)
    return df

def load_file_chunked(path, chunksize=1000000, filter_col=None, filter_values=None, sep=None):
    """
    Load a file in chunks, optionally filtering rows where `filter_col` is in `filter_values`.
    Returns a concatenated DataFrame of valid rows.
    """
    import tqdm
    
    print(f"  Reading {path} in chunks (size={int(chunksize)})...")
    
    # helper: detect sep if not provided
    if sep is None:
        try:
            with open(path, 'r') as f:
                header = f.readline()
                if '\t' in header: sep = '\t'
                elif ',' in header: sep = ','
                else: sep = r'\s+'
        except:
            sep = r'\s+'

    # First, read just the header to get columns and standardize them maps
    # This is tricky because we need to map the filter_col to the raw column name
    # So we read 1 row first
    try:
        header_df = pd.read_csv(path, sep=sep, nrows=1)
    except:
        header_df = pd.read_csv(path, sep=r'\s+', nrows=1) # Fallback
        
    # We want to know which column in the FILE corresponds to our standardized 'filter_col' (e.g. 'SNP')
    # So we run standardize_columns on the header_df to get the mapping
    # But standardize_columns modifies df in place/returns df. We need the mapping dict.
    # Let's slightly refactor or just run it and compare columns.
    
    std_header_df = standardize_columns(header_df.copy())
    
    # Create a map of {RawName: StandardName}
    # This might be hard if standardize_columns doesn't return the map.
    # We can infer it by index if order is preserved, which it matches.
    # Actually, let's just use standardize_columns on the chunks.
    # BUT we need to know which raw column is filter_col to filter EFFICIENTLY? 
    # Pandas read_csv doesn't filter on load easily without an iterator. 
    # We will load chunk -> rename -> filter -> accumulate.
    
    if filter_values is not None:
        filter_set = set(filter_values)
    
    chunks = []
    total_processed = 0
    total_kept = 0
    
    # We use an iterator
    # Note: fallback sep logic for iterator might be fragile if we guessed wrong above.
    # We'll use the one we detected.
    
    try:
        reader = pd.read_csv(path, sep=sep, chunksize=chunksize)
    except:
         reader = pd.read_csv(path, sep=r'\s+', chunksize=chunksize)
         
    # Need to approximate total size for tqdm?
    file_size_mb = os.path.getsize(path) / 1024 / 1024
    
    pbar = tqdm.tqdm(desc="Processing chunks", unit="chunk")
    
    for chunk in reader:
        chunk = standardize_columns(chunk)
        
        # Filter if requested
        if filter_col and filter_values is not None:
            if filter_col in chunk.columns:
                chunk = chunk[chunk[filter_col].isin(filter_set)]
        
        if not chunk.empty:
            chunks.append(chunk)
            total_kept += len(chunk)
            
        total_processed += chunksize
        pbar.update(1)
        
    pbar.close()
    
    if chunks:
        final_df = pd.concat(chunks, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=std_header_df.columns)
        
    print(f"  Finished. Kept {total_kept} rows.")
    return final_df

# Global cache for SNP indices to avoid repeated scans of large BIM files
_SNP_INDEX_CACHE = {}
_SNP_MAP_CACHE = {} # Cache for {snp: index} from BIM reader

def load_bim_with_progress(bim_path):
    """
    Load PLINK .bim file with a progress bar and return a SNP->Index map.
    This is much faster than G.snp.values for huge files because we control the read.
    """
    import tqdm
    if not os.path.exists(bim_path):
        return {}
        
    print(f"  Indexing reference panel: {os.path.basename(bim_path)}... (may take a moment)")
    
    # Estimate size
    size_mb = os.path.getsize(bim_path) / 1024 / 1024
    
    snp_map = {}
    idx_counter = 0
    
    # BIM format: CHROM IF SNP POS A1 A2 (Standard 6 columns)
    # We only need SNP (col 1, 0-indexed) and maybe CHROM if we want to filter earlier.
    # Actually, pandas_plink loads everything into one single dask array.
    # The index in the dask array corresponds to line number in BIM.
    
    # Chunked reading
    chunksize = 1000000
    try:
        # Use C engine for speed
        with pd.read_csv(bim_path, sep=r'\s+', header=None, chunksize=chunksize, 
                         usecols=[0, 1], names=['CHROM', 'SNP'], engine='c', dtype={'CHROM': str, 'SNP': str}) as reader:
            for chunk in tqdm.tqdm(reader, desc="  Reading BIM", unit="chunk"):
                # We need global index.
                # Enumerate local chunk
                snps = chunk['SNP'].values
                # Add to dict
                # This dict creation might be the slow part for 10M SNPs?
                # A dict comprehension is fast.
                new_map = dict(zip(snps, range(idx_counter, idx_counter + len(snps))))
                snp_map.update(new_map)
                idx_counter += len(snps)
                
    except Exception as e:
        print(f"  Error reading BIM: {e}")
        return {}
        
    print(f"  Indexed {len(snp_map)} variants.")
    return snp_map

def get_snp_indices(G, snp_list, chrom=None):
    """
    Get integer indices for a list of SNPs efficiently.
    NOW USES: _SNP_MAP_CACHE populated by main.py calling load_bim_with_progress usually.
    OR builds it on the fly if G is passed but cache empty.
    """
    global _SNP_INDEX_CACHE
    global _SNP_MAP_CACHE
    
    # If using separate BIM loader, we expect _SNP_MAP_CACHE to be populated externally if possible.
    # But if not, we fallback to G.snp.values.
    
    # Check if we have a special map for this G object (by ID or path?)
    # Since G doesn't easily store path, we just use object ID for the xarray style,
    # BUT for the new optimizing plan, the user will likely load the map once.
    # Let's check if there is a 'current_global_map' or similar.
    
    # HEURISTIC: If 'chrom' is provided, we assume the user already subsetted G OR
    # we need to filter the map. The map is global (index in the full BED file).
    # IF G was loaded via read_plink1_bin, it corresponds to the full BIM.
    
    # Let's try to use the cache if it matches the size of G
    g_size = G.shape[0] if G.dims[0]=='variant' else G.shape[1] # Variant dim size
    
    # Use a generic cache key 'GLOBAL' if we assume one reference at a time
    # Or just rebuild if missing.
    
    obj_id = id(G)
    
    if obj_id not in _SNP_INDEX_CACHE:
        # Fallback to old slow method OR use map if available?
        # If the map is in _SNP_MAP_CACHE under a known key...
        # Let's assume the caller put it there? No, that's brittle.
        # Let's use the slow method (G.sel) here ONLY if we haven't pre-loaded.
        
        # NOTE: G.sel is slow.
        # Constructing dict from G.snp.values is faster IF G.snp.values is fast.
        # Optimization:
        if chrom is not None:
             G_subset = G.sel(variant=(G.chrom == str(chrom)))
             mapping = {name: i for i, name in enumerate(G_subset.snp.values)}
             # This mapping provides relative indices for G_subset
             _SNP_INDEX_CACHE[obj_id] = (G_subset, mapping)
        else:
             mapping = {name: i for i, name in enumerate(G.snp.values)}
             _SNP_INDEX_CACHE[obj_id] = (G, mapping)
             
    G_target, cache = _SNP_INDEX_CACHE[obj_id]
    
    found_snps = []
    indices = []
    
    # Intersection logic
    # List comprehension is faster than loop
    # Filter snp_list to those in cache
    # But snp_list is a numpy array usually
    
    # Fast intersection using sets if snp_list is large
    if len(snp_list) > 1000:
        candidates = set(snp_list)
        # Iterate cache? No, cache is big. Iterate candidates.
        # Check existence
        for snp in candidates:
            if snp in cache:
                found_snps.append(snp)
                indices.append(cache[snp])
    else:
        for snp in snp_list:
            if snp in cache: # Dict lookup is O(1)
                found_snps.append(snp)
                indices.append(cache[snp])
            
    return G_target, np.array(indices), np.array(found_snps)


def get_ld_matrix(G, snp_list, chrom=None):
    """
    Extract LD matrix efficiently using index caching.
    """
    G_sub, idx, snps = get_snp_indices(G, snp_list, chrom=chrom)
    if len(idx) <= 1:
        if len(idx) == 1:
            return np.array([[1.0]]), snps
        return None, []
        
    # Trigger read only for selected indices (using isel is much faster)
    # Optimization: Read in chunks if idx is large to avoid RAM spike
    n_variants = len(idx)
    chunk_size = 5000 
    
    genotypes_list = []
    
    import math
    n_chunks = math.ceil(n_variants / chunk_size)
    
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, n_variants)
        sub_idx = idx[start:end]
        
        # Read chunk
        g_chunk = G_sub.isel(variant=sub_idx).values.astype(float)
        
        # Orient to (Samples x Variants) locally
        if G_sub.dims[0] == 'variant':
             g_chunk = g_chunk.T
             
        genotypes_list.append(g_chunk)
        
    if genotypes_list:
        genotypes = np.hstack(genotypes_list) # Stack columns (variants)
    else:
        return None, []
    
    # Simple mean imputation
    means = np.nanmean(genotypes, axis=0)
    inds = np.where(np.isnan(genotypes))
    genotypes[inds] = np.take(means, inds[1])
    # Handle columns with all NaNs
    genotypes = np.nan_to_num(genotypes, nan=0.0)
            
    ld_matrix = np.corrcoef(genotypes, rowvar=False)
    if ld_matrix.ndim == 0:
        ld_matrix = np.array([[1.0]])
        
    return ld_matrix, snps

def get_ld_vector(G, target_snp, snp_list, chrom=None):
    """
    Calculate LD (r) between target_snp and others efficiently.
    """
    G_sub, idx_list, snps = get_snp_indices(G, snp_list, chrom=chrom)
    _, idx_target, _ = get_snp_indices(G, [target_snp], chrom=chrom)
    
    if len(idx_target) == 0 or len(idx_list) == 0:
        return None, []
        
    # Extract target genotype (single variant)
    G_target = G_sub.isel(variant=idx_target).values.flatten().astype(float)
    mu_t = np.nanmean(G_target)
    G_target = np.nan_to_num(G_target, nan=mu_t if not np.isnan(mu_t) else 0)
    
    # Extract list genotypes (chunked)
    n_variants = len(idx_list)
    chunk_size = 5000 
    genotypes_list = []
    import math
    n_chunks = math.ceil(n_variants / chunk_size)
    
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, n_variants)
        sub_idx = idx_list[start:end]
        
        g_chunk = G_sub.isel(variant=sub_idx).values.astype(float)
        if G_sub.dims[0] == 'variant': g_chunk = g_chunk.T
        genotypes_list.append(g_chunk)
        
    if genotypes_list:
        G_list = np.hstack(genotypes_list)
    else:
        return None, []
    
    # Mean impute vector-wise
    means = np.nanmean(G_list, axis=0)
    inds = np.where(np.isnan(G_list))
    G_list[inds] = np.take(means, inds[1])
    G_list = np.nan_to_num(G_list, nan=0.0)
            
    # Pearson r vector calculation
    def corr_vec(X, Y):
        X_m = X - np.mean(X)
        Y_m = Y - np.mean(Y, axis=0)
        denom = (np.sqrt(np.sum(X_m**2)) * np.sqrt(np.sum(Y_m**2, axis=0)))
        denom[denom == 0] = 1.0 # Avoid div by zero
        return np.sum(X_m[:, None] * Y_m, axis=0) / denom
        
    try:
        r_vec = corr_vec(G_target, G_list)
        return r_vec, snps
    except:
        return None, []

def clear_snp_cache():
    """Clear the global SNP index cache to free memory."""
    global _SNP_INDEX_CACHE
    _SNP_INDEX_CACHE.clear()

def read_plink1_bin(*args, **kwargs):
    from pandas_plink import read_plink1_bin as _read
    return _read(*args, **kwargs)
