# RVPRS Pipeline: Mathematical Formulation

## 1. Input Data

We observe summary statistics from $K$ cohorts for a gene with $p$ variants:

- **Score statistics** (per cohort $k$): $U_j^{(k)} = G_j^{(k)\top} Y^{(k)}$, the score for variant $j$
- **Score variance** (per cohort $k$): $V_j^{(k)} = \text{Var}(U_j^{(k)})$
- **LD matrix** (per cohort $k$): $\Lambda^{(k)} = G^{(k)\top} G^{(k)}$, the genotype cross-product (G'G)
- **Allele frequency**: $p_j^{(k)}$ (AF of coded allele in cohort $k$)
- **Sample size**: $n_k$ per cohort

Meta-analytic integration sums across cohorts:

$$U_j = \sum_k U_j^{(k)}, \quad V_j = \sum_k V_j^{(k)}, \quad \Lambda = \sum_k \Lambda^{(k)}, \quad n = \max_j \sum_k n_{k,j}$$

Allele frequency is $N$-weighted: $\hat{p}_j = \sum_k n_{k,j} p_j^{(k)} / \sum_k n_{k,j}$.

**Code:** `preprocessing.py:integrate_summaries()` (lines 260–363)
- Scores summed: `u_sum[idx] += summary_df['Tstat']` (line 310)
- Variances summed: `v_sum[idx] += summary_df['var']` (line 311)
- LD summed: `ld_sum = ld_sum + ld_part` (line 337)
- AF N-weighted: `nw_af_sum[idx] += n * af` (line 321), then `AF = nw_af_sum / n_sum` (line 351)
- $n$: `n = max(N_sum)` across variants (used as `n_scalar` in downstream)
- Markers with zero G'G diagonal dropped (lines 340–344)

**LD symmetry detection:** SAIGE step3 LD files store both triangles $(i,j,v)$ and $(j,i,v)$. The pipeline auto-detects this via `_is_symmetric_csr()` (lines 214–221) and skips symmetrisation to avoid doubling off-diagonal entries. Controlled by `load_ld_coo(..., symmetrize=None)` (lines 224–244). Default: auto-detect.

---

## 2. Preprocessing

### 2.1 Collapsing (MAC-based)

Variants with minor allele count $\text{MAC}_j < T$ (default $T = 10$) are collapsed into a single unit via a collapsing matrix $C$ ($m \times p$, where $m \le p$):

$$U_c = C\, U, \quad \Lambda_c = C\, \Lambda\, C^\top$$

For the collapsed group (row 0 of $C$), the entry is $\sum_j G_j$ over the rare variant set.

**Code:** `preprocessing.py:process_genetic_data()` (lines 479–586)
- MAC computed: `mac = 2 * N * maf` (line 516)
- Collapse indices: `idx_collapse = np.where(mac < mac_threshold)` (line 517)
- Collapse matrix $C$: `get_collapsing_matrix_py()` (line 520, defined lines 401–417)
- Score collapsing: `u_c = collapse_matrix @ u` (line 522)
- G'G collapsing: `collapsed_gtg = collapse_matrix @ ld_matrix @ collapse_matrix.T` (line 528)
- Collapsed MAF: N-weighted `sum(maf * n) / sum(n)` (lines 560–561)

**Binary collapse** (`--binary-collapse`): RareEffect caps collapsed genotypes at 1 (carrier indicator). At the summary-stat level, the collapsed row/column of both $\Phi_c$ and $\Lambda_c$ are rescaled so that the diagonal $\approx \text{MAC}_{\text{total}}$ (number of carriers) rather than the additive sum-of-crossproducts.

**Code:** `preprocessing.py:_apply_binary_collapse_adjustment()` (lines 434–476)
- Ratio: `binary_var_ratio = mac_total / gtg_additive_diag` (line 470)
- Scale factor: `alpha = sqrt(binary_var_ratio)` (line 475)
- Applied to both phi (line 539–541) and gtg (lines 546–548)

### 2.2 Genotype Variance (per-variant N correction)

The genotype variance is computed using per-variant sample sizes $n_j$ for the diagonal, falling back to $n = \max_j n_j$ for off-diagonal normalisation:

$$\text{Var}(G_j) = \frac{\Lambda_{jj}}{n_j} - (2\hat{p}_j)^2$$

**Code:** `preprocessing.py:get_cor_sparse_py()` (lines 366–383)
- Per-variant N passed: `n_per_variant` parameter (line 367)
- Diagonal uses per-variant N: `n_diag = n_per_variant if n_per_variant is not None else n` (line 376)
- Off-diagonal normalised by scalar $n$: `first = (var_div_sparse @ gtg_sparse @ var_div_sparse) / n` (line 381)
- Called with per-variant N from `process_genetic_data()`: line 508

### 2.3 Score Statistic Covariance (phi)

Let $D = \text{diag}(1/\sqrt{\text{Var}(G_j)})$ and define:

$$\Phi_1 = D_U \left(\frac{D\, \Lambda\, D}{n}\right) D_U, \quad \Phi_2 = \frac{2\hat{p}_j}{\sqrt{\text{Var}(G_j)}} \cdot \sqrt{V_j}$$

where $D_U = \text{diag}(\sqrt{V_j})$. After collapsing and mean-centering:

$$\Phi_c = C\, \Phi_1\, C^\top - (C\, \Phi_2)(C\, \Phi_2)^\top$$

This gives $\Phi_c \approx n\, \sigma_Y^2\, \text{Cov}(G_c)$, the score statistic covariance.

**Code:** `preprocessing.py:process_genetic_data()` (lines 512–525)
- `cor_g = cor_s['First']` ← $D \Lambda D / n$ (line 509, from `get_cor_sparse_py`)
- `mean_adj = cor_s['Second_Vec']` ← $2\hat{p}_j / \sqrt{\text{Var}(G_j)}$ (line 510)
- `phi_w1 = dmat @ cor_g @ dmat` ← $\Phi_1 = D_U \cdot (D\Lambda D/n) \cdot D_U$ (line 513)
- `phi_w2 = mean_adj * sd_u` ← $\Phi_2$ (line 514)
- `collapsed_phi = phi_c_w1 - outer(phi_c_w2, phi_c_w2)` ← $\Phi_c$ (line 525)

### 2.4 Genotype Correlation

$$R = D_\Phi^{-1}\, \Phi_c\, D_\Phi^{-1}, \quad D_\Phi = \text{diag}(\sqrt{\text{diag}(\Phi_c)})$$

This recovers $R \approx \text{Cor}(G_c)$, the genotype correlation matrix.

**Code:** `preprocessing.py:covariance_to_correlation()` (lines 386–399)
- `inv_sd = 1 / sqrt(variances)` (line 391)
- `corr = dmat @ cov_sparse @ dmat` (line 393)
- Called at line 580: `collapsed_corr = covariance_to_correlation(collapsed_phi)`

---

## 3. Phenotype Variance Estimation

The profile likelihood (Mode B) requires $Y^\top Y / n = \sigma_Y^2$. When not provided externally, this is estimated from summary statistics:

$$\hat\sigma_Y^2 = \text{median}_j \frac{V_j}{\Lambda_{jj}}$$

**Derivation:** Under the null model, $V_j = \text{Var}(U_j) \approx n \cdot \text{Var}(G_j) \cdot \sigma_Y^2$ and $\Lambda_{jj} \approx n \cdot \text{Var}(G_j)$. The ratio $V_j / \Lambda_{jj} = \sigma_Y^2$ is independent of allele frequency, giving a robust per-variant estimator.

The median is taken over variants with $\Lambda_{jj} > 10$ to exclude ultra-rare variants with unstable estimates. Empirically this gives tight IQR (e.g., $[1.027, 1.028]$ for APOB).

**Code:** `preprocessing.py:estimate_yty_over_n()` (lines 614–633)
- Ratio: `ratios = var_u[mask] / gtg_diag[mask]` (line 629)
- Filter: `mask = (gtg_diag > min_gtg_diag) & (var_u > 0)` (line 625)
- Median: `yty = float(np.median(ratios))` (line 630)
- Auto-invoked in `run_rvprs_pipeline.py:main()` (lines 171–173) when `--genotype-scale covariance` and `--yty-over-n` not provided

**Note:** The HEELS-literature formula $\text{mean}(\hat\beta_j^2 + \text{SE}_j^2 \cdot n)$ is designed for common-variant GWAS with standardised genotypes ($\Lambda_{jj} \approx n$) and does **not** apply to rare variants, where $\Lambda_{jj}$ varies by orders of magnitude across allele frequencies.

---

## 4. Two Estimation Modes

### Mode A: Correlation Scale (HEELS-EM, `--genotype-scale correlation`)

**Inputs to HEELS:**

$$z_m = \frac{U_c}{\sqrt{V_c}} \cdot \sqrt{\frac{n}{m}}, \quad L = R \cdot \frac{n}{m}$$

where $V_c = \text{diag}(\Phi_c)$ are exact collapsed score variances.

**Code:** `preprocessing.py:build_heels_inputs()` (lines 669–683)
- z-score: `z_m = collapsed_score / sqrt(collapsed_var_exact)` (line 672)
- Scaled: `z_m = z_m * sqrt(n / m)` (line 673)
- LD: `ld_corr = ld_corr * n / m` (line 677)

**Model:**

$$z_m = L\, \alpha + \varepsilon, \quad \alpha \sim \mathcal{N}(0,\, \sigma_g\, I_m), \quad \varepsilon \sim \mathcal{N}(0,\, \sigma_e\, L)$$

The eigendecomposition $L = Q\, \text{diag}(d)\, Q^\top$ gives:

$$\hat\alpha = Q\, \text{diag}\!\left(\frac{1}{d_j + \lambda}\right) Q^\top z_m, \quad \lambda = \frac{\sigma_e}{\sigma_g}$$

**Code:** `meta_rvprs_fixed.py:heels_blup_and_trace()` (lines 58–90)
- Projection: `es = ld.eigvec.T @ z_m` (line 73)
- Denominator: `denom = lam + ld.eigval` (line 81)
- BLUP: `beta_hat = ld.eigvec @ (es / denom)` (lines 83–84)

**EM updates** (iterate until convergence):

$$\sigma_g^{(\text{new})} = \frac{\|\hat\alpha\|^2}{m - \text{tr}\!\left((L + \lambda I)^{-1}\right)}, \quad \sigma_e^{(\text{new})} = \frac{Y^\top Y}{n} - \frac{z_m^\top \hat\alpha}{n}$$

When $Y^\top Y / n$ is unknown, renormalization forces $\sigma_g + \sigma_e = 1$.

**Code:** `meta_rvprs_fixed.py:run_heels()` (lines 93–178)
- Lambda: `lam = sigma_e / sigma_g` (line 133)
- sigma_g update: `sigma_g_new = sum(beta_hat**2) / (m - trace_winv)` (lines 136–137)
- sigma_e update: `sigma_e_new = yty_over_n - (z_m @ beta_hat) / n` (line 138)
- Renormalize: `sigma_g_new /= total; sigma_e_new /= total` (lines 142–145)
- Trace: `trace_winv = sum(1 / denom)` (line 88)
- Eigendecomposition: `meta_rvprs_fixed.py:decompose_ld()` (lines 43–55)

**MAF-weighted penalty** (optional, `weights` parameter): replaces the isotropic prior $\alpha \sim \mathcal{N}(0, \sigma_g I)$ with $\alpha \sim \mathcal{N}(0, \sigma_g W)$ where $W = \text{diag}(w_j)$, $w_j = (2p_j(1-p_j) + \varepsilon)^\gamma$. The BLUP becomes:

$$\hat\alpha = (L + \lambda\, W^{-1})^{-1}\, z_m$$

In the eigenbasis of $L$, the penalty projects as $\tilde{w}_j^{-1} = \sum_i Q_{ij}^2 / w_i$ per eigencomponent.

**Code:** `meta_rvprs_fixed.py:heels_blup_and_trace()` (lines 75–79)
- Weight projection: `inv_w_eig = (ld.eigvec.T ** 2) @ (1/weights)` (line 78)
- Weighted denom: `denom = ld.eigval + lam * inv_w_eig` (line 79)
- Passed through `run_heels(..., weights=weights)` (line 106, forwarded at line 134)

**Heritability:**

$$h^2_{\text{corr}} = \frac{\sigma_g}{\sigma_g + \sigma_e}$$

**Code:** `meta_rvprs_fixed.py` line 171: `h2 = sigma_g / (sigma_g + sigma_e)`

**Limitation:** On the correlation scale, all variants are treated as unit-variance features. For mixed-frequency variant sets (e.g., LoF with both rare and common variants), this distorts the relative ranking of effect sizes. The Spearman correlation with RareEffect remains ~0.48 regardless of $\lambda$ or MAF-weighting, because the rank distortion is embedded in the standardised z-scores and LD, not the penalty.

---

### Mode B: Covariance Scale (FaST-LMM MLE, `--genotype-scale covariance`)

**This is the recommended mode for rare variant heritability estimation.**

**Inputs:**

$$U_c \text{ (raw scores)}, \quad \Lambda_c = C\, \Lambda\, C^\top \text{ (collapsed G'G)}, \quad \hat\sigma_Y^2 \text{ (estimated or provided)}$$

**Code:** `preprocessing.py:build_heels_inputs()` (lines 651–667)
- Raw scores: `z_m = summary['collapsed_score']` (line 656) — no z-score transformation
- Raw G'G: `ld = collapsed.collapsed_gtg.toarray()` (line 654) — no n/m scaling

**Model (matching RareEffect):**

$$U_c = \Lambda_c\, \beta + \varepsilon, \quad \beta \sim \mathcal{N}(0,\, \tau\, I_m), \quad \varepsilon \sim \mathcal{N}(0,\, \sigma^2\, \Lambda_c)$$

where $\tau$ is the **per-variant effect variance** on the raw genotype scale (constant for all variants), and $\sigma^2$ is the residual variance.

**Three sufficient statistics from summary data:**

| Quantity | Formula | Code location |
|---|---|---|
| Eigenvalues $S$ of $\Lambda_c$ | $\text{eigh}(\Lambda_c)$ | `meta_rvprs_fixed.py` lines 220–227 |
| Projections $\tilde{Y}_j = S_j^{-1/2} V_j^\top U_c$ | $S^{-1/2} V^\top \cdot \text{score}$ | lines 230–231 |
| Residual $\|Y_\perp\|^2 = n\hat\sigma_Y^2 - \sum_j \tilde{Y}_j^2$ | Requires $\hat\sigma_Y^2$ | lines 234–236 |

**Why this works:** From the SVD $G = U_{\text{svd}}\, S^{1/2}\, V^\top$:

$$U_{\text{svd}}^\top Y = S^{-1/2}\, V^\top\, G^\top Y = S^{-1/2}\, V^\top\, U_c$$

So the projections $\tilde{Y}$ used by RareEffect's `calc_log_lik` are exactly recoverable from score statistics and the G'G eigendecomposition. This is implemented at `meta_rvprs_fixed.py` lines 229–231:
```python
VtU = V.T @ score          # V' G'Y
UtY = VtU / np.sqrt(S)     # S^{-1/2} V' G'Y = U_svd' Y
```

**Profile log-likelihood** (optimised over $\delta = \sigma^2 / \tau$):

$$\ell(\delta) = -\frac{1}{2}\left[n \log(2\pi) + \sum_{j=1}^{k} \log(S_j + \delta) + (n-k)\log\delta + n + n\log\hat\sigma^2(\delta)\right]$$

where

$$\hat\sigma^2(\delta) = \frac{1}{n}\left[\sum_j \frac{\tilde{Y}_j^2}{S_j + \delta} + \frac{\|Y_\perp\|^2}{\delta}\right]$$

**Code:** `meta_rvprs_fixed.py:_fastlmm_profile_loglik()` (lines 181–200)
- `log_lik1 = sum(UtY**2 / (S + delta))` (line 194)
- `log_lik2 = resid_yy / delta` (line 195)
- `sigma_hat = (log_lik1 + log_lik2) / n` (line 196)
- `log_det = sum(log(S + delta)) + (n-k)*log(delta)` (line 199)

Matches RareEffect's `RVPRS_function.R:calc_log_lik()` (lines 188–202) term-by-term.

**Optimisation:** Brent's method over $\delta \in [10^{-10},\, 10^{8}]$, maximising $\ell(\delta)$.

**Code:** `meta_rvprs_fixed.py` lines 241–246:
```python
res = minimize_scalar(lambda d: -_fastlmm_profile_loglik(...), bounds=(1e-10, 1e8), method='bounded')
```

**Variance components:**

$$\hat\sigma^2 = \hat\sigma^2(\hat\delta), \quad \hat\tau = \frac{\hat\sigma^2}{\hat\delta}$$

**Code:** lines 249–252: `sigma_sq = (log_lik1 + log_lik2) / n`, `tau = sigma_sq / opt_delta`

**BLUP:**

$$\hat\beta = ({\Lambda_c + \hat\delta\, I})^{-1}\, U_c = V\, \text{diag}\!\left(\frac{V^\top U_c}{S_j + \hat\delta}\right)$$

**Code:** line 255: `beta_hat = V @ (VtU / (S + opt_delta))`

Matches RareEffect's `RVPRS_function.R:calc_post_beta()` (lines 205–212).

**Heritability:**

$$h^2_{\text{cov}} = \frac{\hat\tau \cdot \text{tr}(\Lambda_c / n)}{\hat\tau \cdot \text{tr}(\Lambda_c / n) + \hat\sigma^2}$$

where $\text{tr}(\Lambda_c / n) = \sum_j \Lambda_{c,jj} / n \approx \sum_j \text{Var}(G_j)$ is the total genotypic variance.

**Code:** lines 258–260:
```python
tr_gtg_over_n = sum(diag(ld_gtg)) / n
h2 = tau * tr_gtg_over_n / (tau * tr_gtg_over_n + sigma_sq)
```

Matches RareEffect's h2 formula in `ComputationEval.R` line 186:
```r
h2_lof_adj <- tau_lof_adj * tr_GtG_lof / (tau_lof_adj * tr_GtG_lof + sigma_sq * n_samples)
```

**Posterior variance:**

$$\text{Var}(\hat\beta_j \mid \text{data}) = \hat\sigma^2 \cdot \left[\sum_\ell \frac{V_{j\ell}^2}{S_\ell + \hat\delta}\right]$$

**Code:** lines 263–264:
```python
posterior_diag = sum((V**2) / (S + opt_delta), axis=1)
posterior_var_diag = sigma_sq * posterior_diag
```

---

## 5. Pipeline Orchestration

**Code:** `run_rvprs_pipeline.py:prepare_and_run_one_annotation()` (lines 34–104)

The function routes to the appropriate estimator based on `genotype_scale`:

```
if genotype_scale == 'covariance':
    yty_val = yty_over_n if yty_over_n is not None else 1.0
    result = run_fastlmm_sumstats(score, ld_gtg, n, yty=yty_val)
else:
    ld = decompose_ld(ld_corr, rank=rank)
    result = run_heels(z_m, ld, m, n, ...)
```

(lines 54–80)

**Auto-estimation of yty:** In `main()` (lines 170–173):
```python
if yty_over_n is None and args.genotype_scale == 'covariance':
    yty_over_n = estimate_yty_over_n(integrated)
```

This is computed once from the full integrated data and shared across all annotations.

---

## 6. Relationship Between the Two Scales

### Eigenvalue relationship

$$d_j^{(\text{corr})} = \frac{n}{m} \cdot r_j \approx \frac{S_j}{m}$$

since $R \approx \Lambda / n$ (for centred genotypes) and $r_j = S_j / n$.

### Lambda relationship

$$\delta_{\text{RE}} = \frac{\sigma^2}{\tau} = m \cdot \lambda_{\text{corr}} = m \cdot \frac{\sigma_e}{\sigma_g}$$

This algebraic identity holds when $\sigma_g = m\,\tau$, but the two models estimate **different** variance components from the same data (different priors, different likelihoods). In practice, $m \cdot \lambda_{\text{corr}}$ can exceed the directly-estimated $\delta$ by orders of magnitude (e.g., 2884× for APOB LoF). **Direct post-hoc conversion of HEELS variance components to covariance-scale h2 does not produce meaningful results.**

### h2 relationship

$$h^2_{\text{corr}} = \frac{\sigma_g}{\sigma_g + \sigma_e} \quad \text{(assumes } \sum_j \text{Var}(G_j) = m \text{, i.e., unit variance)}$$

$$h^2_{\text{cov}} = \frac{\tau \sum_j \text{Var}(G_j)}{\tau \sum_j \text{Var}(G_j) + \sigma^2}$$

These **agree** when $\sum_j \text{Var}(G_j) = m$ (all variants standardised to unit variance).

For rare variants where $\text{Var}(G_j) = 2p_j(1 - p_j) \ll 1$:

$$\sum_j \text{Var}(G_j) \ll m \quad \Longrightarrow \quad h^2_{\text{cov}} < h^2_{\text{corr}}$$

### Prior on raw effect sizes

| | Correlation scale | Covariance scale |
|---|---|---|
| Prior | $\beta_{\text{std},j} \sim \mathcal{N}(0, \sigma_g/m)$ | $\beta_{\text{raw},j} \sim \mathcal{N}(0, \tau)$ |
| Implied raw prior | $\text{Var}(\beta_{\text{raw},j}) = \frac{\sigma_g}{m \cdot \text{Var}(G_j)}$ | $\text{Var}(\beta_{\text{raw},j}) = \tau$ |
| Frequency dependence | Frequency-dependent (rarer = larger) | Frequency-independent |

---

## 7. Effect Size Scaling Between SAIGE and RareEffect

The pipeline's effect sizes differ from RareEffect's by a constant factor:

$$\frac{\hat\beta_{\text{pipeline}}}{\hat\beta_{\text{RareEffect}}} \approx \sqrt{\frac{\sigma_{Y,\text{SAIGE}}^2}{\sigma_{Y,\text{RE}}^2}}$$

This arises because SAIGE score statistics use $\sigma_{Y,\text{SAIGE}}^2 \approx 1.03$ (approximately standardised phenotype), while RareEffect uses GLMM residuals with $\sigma^2_{\text{RE}} = \text{Var}(\text{residuals}) \approx 0.29$. The ratio $\sqrt{1.03 / 0.29} \approx 1.88$ explains the observed scaling.

**Heritability is unaffected** because both $\tau$ and $\sigma^2$ scale with $\sigma_Y^2$, so h2 = $\tau T / (\tau T + \sigma^2)$ is invariant when scores and $Y^\top Y/n$ are internally consistent.

RareEffect's $\sigma^2$: `RareEffect.R` line 56: `sigma_sq <- var(modglmm$residuals)`

---

## 8. Empirical Validation: APOB, UKB WES 470k, f.30780.0.0 (LDL cholesterol)

### Heritability (n ~ 299,265; covariance/MLE with auto-estimated yty = 1.028, binary collapse)

| Annotation | m | Pipeline h2 | RareEffect h2 | Ratio |
|---|---:|---|---|---|
| LoF | 11 | 2.06e-02 | 1.77e-02 | 1.16 |
| Missense | 332 | 9.54e-03 | 1.03e-02 | 0.93 |
| Synonymous | 144 | 1.28e-03 | 6.76e-04 | 1.89 |

### Effect size correlation vs RareEffect

| Annotation | n variants | Pearson | Spearman |
|---|---|---|---|
| LoF | 10 | 0.869 | **0.988** |
| Missense | 331 | **0.997** | **0.994** |
| Synonymous | 143 | 0.889 | 0.953 |

### Effect size correlation: three methods compared

| | RE vs Profile-LL | RE vs HEELS-EM | Profile-LL vs HEELS-EM |
|---|---|---|---|
| LoF Pearson | 0.869 | 0.501 | 0.856 |
| LoF Spearman | **0.988** | 0.479 | 0.491 |
| Missense Pearson | **0.997** | 0.902 | 0.904 |
| Missense Spearman | **0.994** | 0.979 | 0.984 |

Profile-LL (covariance scale) matches RareEffect far better than HEELS-EM (correlation scale), especially for LoF where the mixed-frequency variant set causes rank distortion on the correlation scale.

---

## 9. Implementation Reference

| Component | File | Function | Lines |
|---|---|---|---|
| Meta-analytic integration | `preprocessing.py` | `integrate_summaries()` | 260–363 |
| LD symmetry detection | `preprocessing.py` | `_is_symmetric_csr()` | 214–221 |
| LD loading + auto-detect | `preprocessing.py` | `load_ld_coo()` | 224–244 |
| G'G to correlation | `preprocessing.py` | `get_cor_sparse_py()` | 366–383 |
| Covariance to correlation | `preprocessing.py` | `covariance_to_correlation()` | 386–399 |
| MAC collapsing + binary adj | `preprocessing.py` | `process_genetic_data()` | 479–586 |
| Binary collapse rescaling | `preprocessing.py` | `_apply_binary_collapse_adjustment()` | 439–476 |
| Phenotype variance estimation | `preprocessing.py` | `estimate_yty_over_n()` | 614–633 |
| HEELS input construction | `preprocessing.py` | `build_heels_inputs()` | 636–683 |
| LD eigendecomposition | `meta_rvprs_fixed.py` | `decompose_ld()` | 43–55 |
| HEELS BLUP + trace | `meta_rvprs_fixed.py` | `heels_blup_and_trace()` | 58–90 |
| HEELS-EM iteration | `meta_rvprs_fixed.py` | `run_heels()` | 93–178 |
| Profile log-likelihood | `meta_rvprs_fixed.py` | `_fastlmm_profile_loglik()` | 181–200 |
| FaST-LMM MLE (covariance) | `meta_rvprs_fixed.py` | `run_fastlmm_sumstats()` | 203–279 |
| Pipeline orchestration | `run_rvprs_pipeline.py` | `prepare_and_run_one_annotation()` | 34–104 |
| Auto yty + main loop | `run_rvprs_pipeline.py` | `main()` | 142–219 |
