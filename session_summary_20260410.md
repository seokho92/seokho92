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

### 6. MAF-Weighted Shrinkage (Explored, Not Effective on Correlation Scale)

Attempted two approaches to fix HEELS correlation-scale ranking:
- **Data rescaling** (W^{1/2} L W^{1/2}, W^{1/2} z): crushed rare variant effects to zero
- **Penalty-only** ((L + λW^{-1})^{-1} z): Spearman stuck at 0.479 regardless of γ

Neither fixes the fundamental rank distortion because the z-scores and LD on the correlation scale already encode the wrong signal structure. **Profile-LL on the covariance scale is the correct approach.**

MAF-weighted penalty is available via `--maf-weights GAMMA` for correlation-mode HEELS.

### 7. G'G Matrix Differences: SAIGE LD vs Plink Genotypes

Direct comparison of the collapsed LoF G'G (11×11) between SAIGE step3 LD and RareEffect's plink-derived G'G revealed critical discrepancies in two overlapping indel variants:

#### G'G matrix comparison (APOB LoF)

| Variant | SAIGE G'G diag | RareEffect G'G diag | Match |
|---|---|---|---|
| COLLAPSED (carrier) | 223 | 223 | Yes |
| 2:21006019:CA:C | 54 | 54 | Yes |
| 2:21006087:C:T | 17 | 17 | Yes |
| 2:21006629:AG:A | 16 | 16 | Yes |
| 2:21009304:G:A | 17 | 17 | Yes |
| 2:21010615:G:A | 13 | 13 | Yes |
| 2:21011300:AAC:A | 24 | 24 | Yes |
| 2:21032391:G:A | 41 | 41 | Yes |
| 2:21038086:C:A | 13 | 13 | Yes |
| **2:21043917:GC:G** | **322** | **316** | **No (+6)** |
| **2:21043920:AGCAGCGCG:A** | **323** | **303** | **No (+20)** |

| Off-diagonal | SAIGE | RareEffect | Cauchy-Schwarz bound |
|---|---|---|---|
| **(21043917, 21043920)** | **321** | **299** | 322.5 (SAIGE) / 309.2 (RE) |

These two variants are overlapping indels 3bp apart. SAIGE's dosage-based LD inflates their diagonals and cross-product compared to plink hard-call genotypes. This is the **dominant source** of the remaining h2 discrepancy.

#### Impact on h2 (Profile-LL, same scores, different G'G)

| G'G source | LoF h2 | vs RareEffect |
|---|---|---|
| SAIGE step3 LD | 2.05e-02 | 1.16× |
| **RareEffect plink G'G** | **1.78e-02** | **1.01×** |
| RareEffect reference | 1.77e-02 | 1.00× |

With RareEffect's own G'G, the pipeline matches within 1%.

#### Impact on HEELS (correlation scale)

| G'G source | LoF Pearson | LoF Spearman |
|---|---|---|
| SAIGE step3 LD | 0.501 | 0.479 |
| **RareEffect plink G'G** | **0.889** | **0.915** |

The LoF rank distortion (Spearman 0.48) was **not inherent to the correlation scale** — it was caused by the inflated off-diagonal between the two indels in SAIGE's LD.

### 8. G'G Stabilisation

Implemented `stabilize_gtg()` with three rules to mitigate G'G artifacts from summary statistics:

| Rule | Description | Parameter |
|---|---|---|
| **Rule 2** | Diagonal normalisation: replace G'G_jj with n·2p(1-p) where within 2× of expected (skip collapsed groups) | `--stabilize-gtg` |
| **Rule 3** | Off-diagonal shrinkage: Λ_reg = (1-ρ)Λ + ρ·diag(Λ) | `--gtg-ld-shrinkage RHO` (default 0.01) |
| **Rule 4** | Outlier pair clipping: cap \|Λ_ij / √(Λ_ii Λ_jj)\| at τ | `--gtg-outlier-clip TAU` (default 0.95) |

#### Stabilisation results (APOB, Profile-LL, yty=1.028, binary collapse)

| Annotation | No stabilisation | With stabilisation | RareEffect | Stabilised ratio |
|---|---|---|---|---|
| LoF | 2.06e-02 (1.16×) | **1.60e-02 (0.91×)** | 1.77e-02 | ↓ closer |
| Missense | 9.54e-03 (0.93×) | **9.14e-03 (0.89×)** | 1.03e-02 | ↓ modest |
| Synonymous | 1.28e-03 (1.89×) | 1.22e-03 (1.81×) | 6.76e-04 | ↓ minimal |

For LoF: 1 outlier pair clipped (21043917/21043920), 10/11 diagonals corrected. The h2 drops from 1.16× to 0.91× of RareEffect — a significant improvement from stabilising the two problematic indels.

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
| MAF-weighted penalty | `meta_rvprs_fixed.py` | `--maf-weights GAMMA` for correlation-mode HEELS |
| G'G stabilisation | `preprocessing.py` | `stabilize_gtg()` with diagonal correction, LD shrinkage, outlier clipping |

---

## Final APOB Results

### Profile-LL (covariance mode, yty=1.028, binary collapse, with stabilisation)

| Annotation | Pipeline h2 | RareEffect h2 | Ratio |
|---|---|---|---|
| LoF | 1.60e-02 | 1.77e-02 | 0.91 |
| Missense | 9.14e-03 | 1.03e-02 | 0.89 |
| Synonymous | 1.22e-03 | 6.76e-04 | 1.81 |

### Effect size correlation (Profile-LL, no stabilisation, yty=1.028)

| Annotation | n variants | Pearson | Spearman |
|---|---|---|---|
| LoF | 10 | 0.869 | **0.988** |
| Missense | 331 | **0.997** | **0.994** |
| Synonymous | 143 | 0.889 | 0.953 |

### Three-method comparison (effect size correlation vs RareEffect)

| | RE vs Profile-LL | RE vs HEELS | RE vs HEELS+MAF(γ=1) |
|---|---|---|---|
| LoF Pearson | **0.869** | 0.501 | 0.432 |
| LoF Spearman | **0.988** | 0.479 | 0.479 |
| Missense Pearson | **0.997** | 0.902 | 0.876 |
| Missense Spearman | **0.994** | 0.979 | 0.965 |

---

## Usage

```bash
# Recommended: covariance mode with stabilisation
python -m rvprs_pipeline.run_rvprs_pipeline \
  --gene APOB \
  --groupfile ... \
  --summary part1 part2 ... \
  --marker-info ... \
  --ld ... \
  --outdir results/ \
  --genotype-scale covariance \
  --binary-collapse \
  --stabilize-gtg

# yty_over_n is auto-estimated from Var(U)/G'G_diag
# LD symmetry is auto-detected (no --symmetrize-ld needed)
# G'G stabilisation: diag correction + 1% LD shrinkage + 0.95 outlier clip

# Correlation mode with MAF weights (alternative)
python -m rvprs_pipeline.run_rvprs_pipeline \
  --gene APOB ... \
  --genotype-scale correlation \
  --maf-weights 1.0
```

---

## Files Modified

- `rvprs_pipeline/preprocessing.py` — CollapsedInputs.collapsed_gtg, estimate_yty_over_n(), stabilize_gtg(), auto-detect symmetry, per-variant N, binary collapse for G'G, build_heels_inputs covariance mode + stabilisation params
- `rvprs_pipeline/meta_rvprs_fixed.py` — run_fastlmm_sumstats(), MAF-weighted penalty in heels_blup_and_trace()
- `rvprs_pipeline/run_rvprs_pipeline.py` — --genotype-scale, --maf-weights, --stabilize-gtg, --gtg-ld-shrinkage, --gtg-outlier-clip, auto yty estimation
- `rvprs_pipeline/pipeline_mathematics.md` — full mathematical documentation with code cross-references
