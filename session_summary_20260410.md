# RVPRS Pipeline Session Summary — 2026-04-09/10

## Objective

Investigate and resolve heritability estimation discrepancies between the RVPRS summary-statistics pipeline and RareEffect (individual-level FaST-LMM), using APOB / LDL cholesterol (f.30780.0.0, UKB WES 470k, n ≈ 299,265) as the test case.

---

## Key Findings

### 1. Genotype Scale Mismatch (Root Cause of h2 Discrepancy)

The pipeline operated on the **correlation scale** (standardised genotypes, HEELS-EM), while RareEffect operates on the **covariance scale** (raw G'G, FaST-LMM MLE). These impose different priors:

| | Correlation (HEELS) | Covariance (RareEffect) |
|---|---|---|
| Prior | Var(β_raw) ∝ 1/Var(G_j) | Var(β_raw) = τ (constant) |
| h2 formula | σ_g / (σ_g + σ_e) | τ · tr(G'G/n) / (τ · tr(G'G/n) + σ²) |
| Rare variant treatment | Inflated prior (upweighted) | Equal prior (frequency-independent) |

### 2. HEELS-EM Cannot Solve on the Covariance Scale

Directly feeding raw G'G into HEELS-EM caused σ_g to collapse to zero for missense/synonymous. The EM fails because `U'β̂/n` is negligible (~100/299k), pinning σ_e ≈ 1 while σ_g spirals down.

**Solution**: Implemented **FaST-LMM profile-likelihood MLE** (Brent optimisation over δ = σ²/τ), matching RareEffect's `fast_lmm()` exactly. This works from summary statistics using the identity:

```
UtY = S^{-1/2} V' (G'Y)     ← recovered from scores and G'G eigendecomposition
resid_yy = n · yty - Σ UtY²  ← requires Y'Y/n
```

### 3. LD Symmetrisation Bug (`--symmetrize-ld`)

SAIGE step3 LD files store **both triangles** `(i,j,v)` and `(j,i,v)`. The `--symmetrize-ld` flag applied `ld + ld.T - diag(ld)`, which **doubled all off-diagonal entries**. This produced a non-PSD G'G with a negative eigenvalue of -319.5 (30% of trace for LoF).

**Fix**: Auto-detect symmetry (`_is_symmetric_csr()`). Default changed from `symmetrize=False` to `symmetrize=None` (auto-detect, skip if already symmetric).

### 4. Y'Y/n Estimation from Summary Statistics

The profile likelihood requires Y'Y/n, unavailable from standard summary stats. Two approaches evaluated:

| Method | Formula | Result | Applicable to |
|---|---|---|---|
| HEELS notebook | mean(β² + SE²·n) | 32,097 (wrong) | Common variants only |
| **Var(U)/G'G_diag** | **median(Var(U_j) / G'G_jj)** | **1.028** | **Any allele frequency** |

The Var(U)/G'G_diag estimator works because Var(U_j) = n·Var(G_j)·σ_Y² and G'G_jj = n·Var(G_j), so Var(G_j) cancels. Consistent across all annotations (IQR: [1.027, 1.028]).

Implemented as `estimate_yty_over_n()` in preprocessing.py, auto-invoked in covariance mode.

### 5. Effect Size Comparison

The ~1.88× scaling between pipeline and RareEffect effect sizes is explained by phenotype variance: SAIGE uses σ_Y² ≈ 1.028, RareEffect uses GLMM residuals σ² ≈ 0.292. The ratio sqrt(1.028/0.292) = 1.876 matches perfectly.

### 6. MAF-Weighted Shrinkage (Explored, Not Effective)

Attempted two approaches to fix HEELS correlation-scale ranking:
- **Data rescaling** (W^{1/2} L W^{1/2}, W^{1/2} z): crushed rare variant effects to zero
- **Penalty-only** ((L + λW^{-1})^{-1} z): Spearman stuck at 0.479 regardless of γ

Neither fixes the fundamental rank distortion because the z-scores and LD on the correlation scale already encode the wrong signal structure. **Profile-LL on the covariance scale is the correct approach.**

---

## Bugs Fixed

| Bug | Location | Fix |
|---|---|---|
| LD symmetrisation doubling off-diagonals | `load_ld_coo()` | Auto-detect symmetry, skip if already symmetric |
| Per-variant N in variance computation | `get_cor_sparse_py()` | Accept `n_per_variant` array for diagonal |
| Collapsed variant MAF | `process_genetic_data()` | N-weighted AF instead of naive mean |
| Incorrect h2_rareeffect_scale conversion | `run_rvprs_pipeline.py` | Removed; replaced by Profile-LL |
| Binary collapse not applied to G'G | `process_genetic_data()` | Apply `_apply_binary_collapse_adjustment` to collapsed_gtg |

---

## New Features

| Feature | Location | Description |
|---|---|---|
| `--genotype-scale covariance` | `run_rvprs_pipeline.py` | FaST-LMM MLE on raw G'G scale, aligned with RareEffect |
| `run_fastlmm_sumstats()` | `meta_rvprs_fixed.py` | Profile-likelihood MLE from scores + G'G |
| `estimate_yty_over_n()` | `preprocessing.py` | Auto-estimate σ_Y² from Var(U)/G'G_diag |
| Auto-detect LD symmetry | `preprocessing.py` | `_is_symmetric_csr()`, default `symmetrize=None` |
| MAF-weighted penalty | `meta_rvprs_fixed.py` | `weights` parameter in `heels_blup_and_trace()` / `run_heels()` |

---

## Final APOB Results (Profile-LL, yty=1.028, binary collapse)

| Annotation | Pipeline h2 | RareEffect h2 | Ratio | Pearson r(β) | Spearman r(β) |
|---|---|---|---|---|---|
| LoF | 2.06e-02 | 1.77e-02 | 1.16 | 0.869 | **0.988** |
| Missense | 9.54e-03 | 1.03e-02 | 0.93 | **0.997** | **0.994** |
| Synonymous | 1.28e-03 | 6.76e-04 | 1.89 | 0.889 | 0.953 |

---

## Usage

```bash
# RareEffect-aligned mode (recommended for rare variant h2)
python -m rvprs_pipeline.run_rvprs_pipeline \
  --gene APOB \
  --groupfile ... \
  --summary part1 part2 ... \
  --marker-info ... \
  --ld ... \
  --outdir results/ \
  --genotype-scale covariance \
  --binary-collapse

# yty_over_n is auto-estimated from Var(U)/G'G_diag
# LD symmetry is auto-detected (no --symmetrize-ld needed)
```

---

## Files Modified

- `rvprs_pipeline/preprocessing.py` — CollapsedInputs.collapsed_gtg, estimate_yty_over_n(), auto-detect symmetry, per-variant N, binary collapse for G'G, build_heels_inputs covariance mode
- `rvprs_pipeline/meta_rvprs_fixed.py` — run_fastlmm_sumstats(), MAF-weighted penalty in heels_blup_and_trace()
- `rvprs_pipeline/run_rvprs_pipeline.py` — --genotype-scale, auto yty estimation, genotype_scale passthrough
- `rvprs_pipeline/pipeline_mathematics.md` — full mathematical documentation (new file)
