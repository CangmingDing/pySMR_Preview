import numpy as np
from scipy import stats

def smr_test(b_zx, se_zx, b_zy, se_zy):
    """
    Perform SMR test.
    Returns: b_xy, se_xy, p_smr
    """
    b_xy = b_zy / b_zx
    # Delta method for SE
    se_xy = np.sqrt((se_zy**2 * b_zx**2 + se_zx**2 * b_zy**2) / (b_zx**4))
    
    z_smr = b_xy / se_xy
    p_smr = stats.chi2.sf(z_smr**2, 1)
    
    return b_xy, se_xy, p_smr

def pchisqsum(x, lambda_vals):
    """
    P-value for a weighted sum of chi-squared variables sum(lambda_i * chi2_1).
    Uses Satterthwaite approximation.
    """
    if len(lambda_vals) == 0:
        return 1.0
    
    sum_lambda = np.sum(lambda_vals)
    sum_lambda_sq = np.sum(lambda_vals**2)
    
    if sum_lambda <= 0:
        return 1.0
        
    a = sum_lambda_sq / sum_lambda
    df = (sum_lambda**2) / sum_lambda_sq
    
    return stats.chi2.sf(x / a, df)

def heidi_test(b_zx, se_zx, b_zy, se_zy, ld_matrix):
    """
    Perform HEIDI test for heterogeneity.
    b_zx, se_zx, b_zy, se_zy: vectors for multiple SNPs
    ld_matrix: LD matrix for these SNPs
    """
    nsnp = len(b_zx)
    if nsnp <= 1:
        return np.nan, 0
    
    b_xy = b_zy / b_zx
    z_zx = b_zx / se_zx
    
    # 1. Identify the 'top' SNP (instrument)
    # SMR v1.4.0 uses max abs z_zx if not specified
    maxid = np.argmax(np.abs(z_zx))
    
    # 2. Estimate covariance of b_xy
    # Formula from est_cov_bxy:
    # cov(b_xy_i, b_xy_j) = LD_ij * (se_zy_i * se_zy_j / (b_zx_i * b_zx_j) + b_xy_i * b_xy_j / (z_zx_i * z_zx_j)) 
    #                       - (b_xy_i * b_xy_j) / (z_zx_i * z_zx_j * z_zx_i * z_zx_j) -- wait, the C++ code is a bit different.
    
    # From C++:
    # bxytbxy = _bxy * _bxy.T
    # zsxztzsxz = _zsxz * _zsxz.T
    # covbxy = _LD_heidi * ((_seyz * _seyz.T) / (_bxz * _bxz.T) + bxytbxy / zsxztzsxz) - bxytbxy / (zsxztzsxz * zsxztzsxz)
    
    bxytbxy = np.outer(b_xy, b_xy)
    zsxztzsxz = np.outer(z_zx, z_zx)
    seyz_outer = np.outer(se_zy, se_zy)
    bxz_outer = np.outer(b_zx, b_zx)
    
    cov_bxy = ld_matrix * (seyz_outer / bxz_outer + bxytbxy / zsxztzsxz) - bxytbxy / (zsxztzsxz**2)
    
    # 3. Calculate deviations from the top SNP
    # d_i = b_xy_top - b_xy_i
    # We want the distribution of D' V^-1 D
    
    # Deviations
    indices = np.delete(np.arange(nsnp), maxid)
    d = b_xy[maxid] - b_xy[indices]
    
    # Variance of deviations:
    # Var(b_xy_top - b_xy_i) = Var(b_xy_top) + Var(b_xy_i) - 2*Cov(b_xy_top, b_xy_i)
    # The full covariance matrix of d is:
    # V_ij = Cov(b_xy_top - b_xy_i, b_xy_top - b_xy_j)
    #      = Var(b_xy_top) + Cov(b_xy_i, b_xy_j) - Cov(b_xy_top, b_xy_i) - Cov(b_xy_top, b_xy_j)
    
    v_top = cov_bxy[maxid, maxid]
    cov_indices = cov_bxy[np.ix_(indices, indices)]
    cov_top_indices = cov_bxy[maxid, indices]
    
    # Matrix V
    V = v_top + cov_indices - np.outer(np.ones(nsnp-1), cov_top_indices) - np.outer(cov_top_indices, np.ones(nsnp-1))
    
    # Add small ridge for stability
    V += np.eye(nsnp-1) * 1e-8
    
    # Chi-square sum statistic
    # sum_chisq = sum(d_i^2 / var(d_i)) -- this is not quite right if they are correlated.
    # SMR actually computes eigenvalues of the correlation matrix of deviations?
    # Let's re-read:
    # sumChisq_dev += dev[i]*dev[i] / tmp3[i] where tmp3 is diagonal of V
    # Then it calculates eigenvalues of correlation matrix of V
    
    diag_V = np.diag(V)
    sum_chisq_dev = np.sum(d**2 / diag_V)
    
    # Correlation matrix of deviations
    V_corr = V / np.sqrt(np.outer(diag_V, diag_V))
    
    # Eigenvalues
    lambdas = np.linalg.eigvalsh(V_corr)
    
    # P-value
    p_heidi = pchisqsum(sum_chisq_dev, lambdas)
    
    return p_heidi, nsnp
