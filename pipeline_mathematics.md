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

**LD symmetry detection:** SAIGE step3 LD files store both triangles $(i,j,v)$ and $(j,i,v)$. The pipeline auto-detects this and skips symmetrisation to avoid doubling off-diagonal entries.

---

## 2. Preprocessing

### 2.1 Collapsing (MAC-based)

Variants with minor allele count $\text{MAC}_j < T$ (default $T = 10$) are collapsed into a single unit via a collapsing matrix $C$ ($m \times p$, where $m \le p$):

$$U_c = C\, U, \quad \Lambda_c = C\, \Lambda\, C^\top$$

For the collapsed group (row 0 of $C$), the entry is $\sum_j G_j$ over the rare variant set.

**Binary collapse** (`--binary-collapse`): RareEffect caps collapsed genotypes at 1 (carrier indicator). At the summary-stat level, the collapsed row/column of $\Lambda_c$ is rescaled so that $\Lambda_{c,00} \approx \text{MAC}_{\text{total}}$ (number of carriers) rather than the additive sum-of-crossproducts.

### 2.2 Genotype Variance (per-variant N correction)

The genotype variance is computed using per-variant sample sizes $n_j$ for the diagonal, falling back to $n = \max_j n_j$ for off-diagonal normalisation:

$$\text{Var}(G_j) = \frac{\Lambda_{jj}}{n_j} - (2\hat{p}_j)^2$$

### 2.3 Score Statistic Covariance (phi)

Let $D = \text{diag}(1/\sqrt{\text{Var}(G_j)})$ and define:

$$\Phi_1 = D_U \left(\frac{D\, \Lambda\, D}{n}\right) D_U, \quad \Phi_2 = \frac{2\hat{p}_j}{\sqrt{\text{Var}(G_j)}} \cdot \sqrt{V_j}$$

where $D_U = \text{diag}(\sqrt{V_j})$. After collapsing and mean-centering:

$$\Phi_c = C\, \Phi_1\, C^\top - (C\, \Phi_2)(C\, \Phi_2)^\top$$

This gives $\Phi_c \approx n\, \sigma_Y^2\, \text{Cov}(G_c)$, the score statistic covariance.

### 2.4 Genotype Correlation

$$R = D_\Phi^{-1}\, \Phi_c\, D_\Phi^{-1}, \quad D_\Phi = \text{diag}(\sqrt{\text{diag}(\Phi_c)})$$

This recovers $R \approx \text{Cor}(G_c)$, the genotype correlation matrix.

---

## 3. Phenotype Variance Estimation

The profile likelihood (Mode B) requires $Y^\top Y / n = \sigma_Y^2$. When not provided externally, this is estimated from summary statistics:

$$\hat\sigma_Y^2 = \text{median}_j \frac{V_j}{\Lambda_{jj}}$$

**Derivation:** Under the null model, $V_j = \text{Var}(U_j) \approx n \cdot \text{Var}(G_j) \cdot \sigma_Y^2$ and $\Lambda_{jj} \approx n \cdot \text{Var}(G_j)$. The ratio $V_j / \Lambda_{jj} = \sigma_Y^2$ is independent of allele frequency, giving a robust per-variant estimator.

The median is taken over variants with $\Lambda_{jj} > 10$ to exclude ultra-rare variants with unstable estimates. Empirically this gives tight IQR (e.g., $[1.027, 1.028]$ for APOB).

**Note:** The HEELS-literature formula $\text{mean}(\hat\beta_j^2 + \text{SE}_j^2 \cdot n)$ is designed for common-variant GWAS with standardised genotypes and does **not** apply to rare variants, where $V_j$ varies by orders of magnitude across allele frequencies.

---

## 4. Two Estimation Modes

### Mode A: Correlation Scale (HEELS-EM, `--genotype-scale correlation`)

**Inputs to HEELS:**

$$z_m = \frac{U_c}{\sqrt{V_c}} \cdot \sqrt{\frac{n}{m}}, \quad L = R \cdot \frac{n}{m}$$

where $V_c = \text{diag}(\Phi_c)$ are exact collapsed score variances.

**Model:**

$$z_m = L\, \alpha + \varepsilon, \quad \alpha \sim \mathcal{N}(0,\, \sigma_g\, I_m), \quad \varepsilon \sim \mathcal{N}(0,\, \sigma_e\, L)$$

The eigendecomposition $L = Q\, \text{diag}(d)\, Q^\top$ gives:

$$\hat\alpha = Q\, \text{diag}\!\left(\frac{1}{d_j + \lambda}\right) Q^\top z_m, \quad \lambda = \frac{\sigma_e}{\sigma_g}$$

**EM updates** (iterate until convergence):

$$\sigma_g^{(\text{new})} = \frac{\|\hat\alpha\|^2}{m - \text{tr}\!\left((L + \lambda I)^{-1}\right)}, \quad \sigma_e^{(\text{new})} = \frac{Y^\top Y}{n} - \frac{z_m^\top \hat\alpha}{n}$$

When $Y^\top Y / n$ is unknown, renormalization forces $\sigma_g + \sigma_e = 1$.

**MAF-weighted penalty** (optional, `weights` parameter): replaces the isotropic prior $\alpha \sim \mathcal{N}(0, \sigma_g I)$ with $\alpha \sim \mathcal{N}(0, \sigma_g W)$ where $W = \text{diag}(w_j)$, $w_j = (2p_j(1-p_j) + \varepsilon)^\gamma$. The BLUP becomes:

$$\hat\alpha = (L + \lambda\, W^{-1})^{-1}\, z_m$$

In the eigenbasis of $L$, the penalty projects as $\tilde{w}_j^{-1} = \sum_i Q_{ij}^2 / w_i$ per eigencomponent.

**Heritability:**

$$h^2_{\text{corr}} = \frac{\sigma_g}{\sigma_g + \sigma_e}$$

This is h2 under a **frequency-dependent** prior on raw effects:
$\text{Var}(\beta_{\text{raw},j}) = \sigma_g / (m \cdot \text{Var}(G_j))$.

**Limitation:** On the correlation scale, all variants are treated as unit-variance features. For mixed-frequency variant sets (e.g., LoF with both rare and common variants), this distorts the relative ranking of effect sizes. The Spearman correlation with RareEffect remains ~0.48 regardless of $\lambda$ or MAF-weighting, because the rank distortion is embedded in the standardised z-scores and LD, not the penalty.

---

### Mode B: Covariance Scale (FaST-LMM MLE, `--genotype-scale covariance`)

**This is the recommended mode for rare variant heritability estimation.**

**Inputs:**

$$U_c \text{ (raw scores)}, \quad \Lambda_c = C\, \Lambda\, C^\top \text{ (collapsed G'G)}, \quad \hat\sigma_Y^2 \text{ (estimated or provided)}$$

**Model (matching RareEffect):**

$$U_c = \Lambda_c\, \beta + \varepsilon, \quad \beta \sim \mathcal{N}(0,\, \tau\, I_m), \quad \varepsilon \sim \mathcal{N}(0,\, \sigma^2\, \Lambda_c)$$

where $\tau$ is the **per-variant effect variance** on the raw genotype scale (constant for all variants), and $\sigma^2$ is the residual variance.

**Three sufficient statistics from summary data:**

| Quantity | Formula | Source |
|---|---|---|
| Eigenvalues $S$ of $\Lambda_c$ | $\text{eigh}(\Lambda_c)$ | G'G (LD files) |
| Projections $\tilde{Y}_j = S_j^{-1/2} V_j^\top U_c$ | $S^{-1/2} V^\top \cdot \text{score}$ | Score stats + G'G |
| Residual $\|Y_\perp\|^2 = n\hat\sigma_Y^2 - \sum_j \tilde{Y}_j^2$ | Requires $\hat\sigma_Y^2$ | Estimated from $V_j / \Lambda_{jj}$ |

**Why this works:** From the SVD $G = U_{\text{svd}}\, S^{1/2}\, V^\top$:

$$U_{\text{svd}}^\top Y = S^{-1/2}\, V^\top\, G^\top Y = S^{-1/2}\, V^\top\, U_c$$

So the projections $\tilde{Y}$ used by RareEffect's `calc_log_lik` are exactly recoverable from score statistics and the G'G eigendecomposition.

**Profile log-likelihood** (optimised over $\delta = \sigma^2 / \tau$):

$$\ell(\delta) = -\frac{1}{2}\left[n \log(2\pi) + \sum_{j=1}^{k} \log(S_j + \delta) + (n-k)\log\delta + n + n\log\hat\sigma^2(\delta)\right]$$

where

$$\hat\sigma^2(\delta) = \frac{1}{n}\left[\sum_j \frac{\tilde{Y}_j^2}{S_j + \delta} + \frac{\|Y_\perp\|^2}{\delta}\right]$$

**Optimisation:** Brent's method over $\delta \in [10^{-10},\, 10^{8}]$, maximising $\ell(\delta)$.

**Variance components:**

$$\hat\sigma^2 = \hat\sigma^2(\hat\delta), \quad \hat\tau = \frac{\hat\sigma^2}{\hat\delta}$$

**BLUP:**

$$\hat\beta = ({\Lambda_c + \hat\delta\, I})^{-1}\, U_c = V\, \text{diag}\!\left(\frac{V^\top U_c}{S_j + \hat\delta}\right)$$

**Heritability:**

$$h^2_{\text{cov}} = \frac{\hat\tau \cdot \text{tr}(\Lambda_c / n)}{\hat\tau \cdot \text{tr}(\Lambda_c / n) + \hat\sigma^2}$$

where $\text{tr}(\Lambda_c / n) = \sum_j \Lambda_{c,jj} / n \approx \sum_j \text{Var}(G_j)$ is the total genotypic variance.

**Interpretation:** Constant per-variant effect prior (frequency-independent). Rarer variants do NOT get inflated priors. This matches RareEffect's FaST-LMM model exactly.

---

## 5. Relationship Between the Two Scales

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

## 6. Posterior Variance

In both modes, the posterior variance of $\hat\beta_j$ is:

$$\text{Var}(\hat\beta_j \mid \text{data}) = \sigma_{\text{noise}} \cdot \left[\sum_\ell \frac{V_{j\ell}^2}{S_\ell + \delta}\right]$$

where $V$ are eigenvectors and $\sigma_{\text{noise}}$ is the relevant noise variance ($\sigma_e$ in correlation mode, $\hat\sigma^2$ in covariance mode).

---

## 7. Effect Size Scaling Between SAIGE and RareEffect

The pipeline's effect sizes differ from RareEffect's by a constant factor:

$$\frac{\hat\beta_{\text{pipeline}}}{\hat\beta_{\text{RareEffect}}} \approx \sqrt{\frac{\sigma_{Y,\text{SAIGE}}^2}{\sigma_{Y,\text{RE}}^2}}$$

This arises because SAIGE score statistics use $\sigma_{Y,\text{SAIGE}}^2 \approx 1.03$ (approximately standardised phenotype), while RareEffect uses GLMM residuals with $\sigma^2_{\text{RE}} = \text{Var}(\text{residuals}) \approx 0.29$. The ratio $\sqrt{1.03 / 0.29} \approx 1.88$ explains the observed scaling.

**Heritability is unaffected** because both $\tau$ and $\sigma^2$ scale with $\sigma_Y^2$, so h2 = $\tau T / (\tau T + \sigma^2)$ is invariant when scores and $Y^\top Y/n$ are internally consistent.

---

## 8. Empirical Validation: APOB, UKB WES 470k, f.30780.0.0 (LDL cholesterol)

### Heritability (n ~ 299,265; covariance/MLE with yty = 1.028, binary collapse)

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

| Component | File | Function |
|---|---|---|
| Meta-analytic integration | `preprocessing.py` | `integrate_summaries()` |
| LD symmetry detection | `preprocessing.py` | `_is_symmetric_csr()` |
| MAC collapsing + binary adjustment | `preprocessing.py` | `process_genetic_data()` |
| G'G to correlation | `preprocessing.py` | `get_cor_sparse_py()` |
| Phenotype variance estimation | `preprocessing.py` | `estimate_yty_over_n()` |
| HEELS input construction | `preprocessing.py` | `build_heels_inputs()` |
| HEELS-EM (correlation mode) | `meta_rvprs_fixed.py` | `run_heels()` |
| FaST-LMM MLE (covariance mode) | `meta_rvprs_fixed.py` | `run_fastlmm_sumstats()` |
| MAF-weighted shrinkage | `meta_rvprs_fixed.py` | `heels_blup_and_trace(weights=...)` |
| Pipeline orchestration | `run_rvprs_pipeline.py` | `prepare_and_run_one_annotation()` |
