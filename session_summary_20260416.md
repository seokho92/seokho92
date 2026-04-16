# RVPRS Pipeline Session Summary — 2026-04-16

## Objective

Build a comprehensive simulation study comparing RVPRS (summary-statistics-based rare-variant PRS) against multiple comparable methods across a systematic scenario grid, following the design principles from `simulation_pipeline.md`.

---

## Background

The existing simulation (`simulate_rvprs_study.py`) compared only three approaches: HEELS-exact (individual-level pseudodata), HEELS-reconstructed (from summary statistics), and Ridge baseline. This was insufficient for a publication-quality evaluation because:

1. No comparison against standard rare-variant PRS alternatives (Burden, IVW, SKAT-O weights)
2. No sparse penalization baselines (LASSO, Elastic Net)
3. No out-of-sample prediction evaluation (all metrics were in-sample)
4. No stress tests (LD mismatch, varying architectures)
5. No systematic scenario grid

The `simulation_pipeline.md` research report surveyed simulation designs from COJO, LDSC, SMR, MTAG, CAUSE, and REMETA, recommending a **two-layer scenario design** (baseline + stress tests) with multiple evaluation metrics.

---

## What Was Built

### New File: `simulation_study.py`

A self-contained simulation study framework (~550 lines) that:

- **Imports and reuses** all data generation functions from `simulate_rvprs_study.py` (genotype simulation, phenotype generation, summary statistic computation, cohort splitting)
- **Imports** HEELS/FaST-LMM estimators from `meta_rvprs_fixed.py`
- **Adds no duplicate logic** — all simulation infrastructure is shared

### 9 Methods Compared

| Method | Data Required | h2 Estimation | Posterior Intervals |
|--------|--------------|:---:|:---:|
| **Ridge** | Individual-level | - | - |
| **LASSO** (5-fold CV) | Individual-level | - | - |
| **Elastic Net** (CV, l1_ratio=[0.1,0.5,0.9]) | Individual-level | - | - |
| **OLS** | Individual-level (n>p only) | - | - |
| **RVPRS-HEELS** | Summary statistics | Yes | Yes |
| **RVPRS-FaST-LMM** | Summary statistics | Yes | Yes |
| **Burden PRS** | Summary statistics | - | - |
| **IVW** (marginal) | Summary statistics | - | - |
| **SKAT-O weights** | Summary statistics | - | - |

Individual-level methods serve as oracle baselines — they use the genotype matrix directly, which summary-statistics methods cannot access in practice. The comparison shows how much information is lost by working from summaries only.

#### Method Details

**RVPRS-HEELS** (our method, correlation scale): Transforms score statistics to z-scores, extracts LD correlation from reference panel, runs HEELS-EM to jointly estimate variance components (h2) and BLUP effects. Back-transforms to original scale.

**RVPRS-FaST-LMM** (our method, covariance scale): Uses score vector and G'G matrix directly, estimates variance components via profile likelihood (Brent optimization over delta = sigma2/tau), computes BLUP. Matches RareEffect's `fast_lmm()`.

**Burden PRS**: Per annotation group (LoF/missense/synonymous), estimates a single group-level effect = sum(U_j) / sum(G'G_jj). All variants in a group share this effect. No per-variant differentiation.

**IVW (marginal)**: Per-variant marginal OLS estimate = U_j / (G'G)_jj. No shrinkage, no LD correction. Simplest possible per-variant estimator.

**SKAT-O weights**: Like IVW but with Beta(MAF; 1, 25) weights that upweight rarer variants, inspired by the SKAT-O kernel weighting scheme.

### 7-Axis Scenario Grid

| Axis | Values | Rationale |
|------|--------|-----------|
| Sample size N | 10k, 50k, 200k | Power scaling |
| Num variants p | 40, 80, 200 | Gene size variation |
| Architecture | sparse, polygenic, mixed | Different genetic models |
| Heritability h2 | 0.01, 0.05, 0.15 | Signal strength |
| Trait type | quantitative, binary (prev=5%) | Phenotype model |
| LD mismatch | matched, mismatched | Reference panel quality stress test |
| Num cohorts | 1, 6 | Meta-analysis structure |

Full grid = 3 x 3 x 3 x 3 x 2 x 2 x 2 = **648 scenarios** x reps.

#### Architecture Presets

The three architectures modify causal variant probability and effect magnitude:

| Architecture | LoF causal% | Missense causal% | Syn causal% | Effect scale |
|-------------|:-----------:|:-----------------:|:-----------:|:------------:|
| **Sparse** | 15% | 3% | 0% | Large (1.5x LoF) |
| **Polygenic** | 80% | 60% | 30% | Small (0.15x all) |
| **Mixed** | 30% | 10% | 2% | Moderate (default) |

#### LD Mismatch Implementation

When `ld_mismatch = "mismatched"`:
- A separate reference panel of size n/5 is generated with different LD structure (rho=0.30 vs generation rho=0.15)
- Summary-stat methods use the reference panel's LD correlation, but the score vector U from the actual training data
- Individual-level methods are unaffected (they don't use the reference panel)
- This simulates real-world scenarios where the LD reference panel comes from a different/smaller cohort

### 8 Evaluation Metrics

| Metric | Formula | What it measures |
|--------|---------|-----------------|
| RMSE | $\sqrt{\text{mean}(\hat\beta - \beta_{\text{true}})^2}$ | Effect estimation accuracy |
| Pearson corr | $\text{cor}(\beta_{\text{true}}, \hat\beta)$ | Linear agreement |
| Spearman corr | $\rho_s(\beta_{\text{true}}, \hat\beta)$ | Rank agreement |
| Bias | $\text{mean}(\hat\beta - \beta_{\text{true}})$ | Systematic over/under-estimation |
| h2 bias | $\hat{h}^2 - h^2_{\text{truth}}$ | Variance component calibration |
| Prediction R^2 | $1 - \text{SS}_{\text{res}} / \text{SS}_{\text{tot}}$ on held-out test set | Out-of-sample prediction |
| 95% CI coverage | Fraction where $|\beta_{\text{true}} - \hat\beta| < 1.96 \cdot \text{posterior\_sd}$ | Interval calibration |
| Wall time | seconds | Computational cost |

**Key design choice**: Train/test split (80/20) applied *before* any summary statistic computation. The test set is used only for prediction R^2, ensuring honest out-of-sample evaluation.

### Output Structure

```
outdir/
  study_config.json           # full config for reproducibility
  replicate_results.tsv       # one row per (scenario, rep, method) — raw data
  aggregated_results.tsv      # one row per (scenario, method) — mean + SD
  <scenario_id>/
    rep_0000_variants.tsv     # per-variant betas across all methods (rep 0)
```

---

## Early Results (Smoke Tests)

### Quantitative trait, h2=0.15, matched LD, N=5000, p=40

| Method | RMSE | Pearson corr | Pred R^2 |
|--------|:----:|:------------:|:--------:|
| LASSO | 0.410 | 0.642 | 0.103 |
| Elastic Net | 0.410 | 0.642 | 0.103 |
| Ridge | 0.439 | 0.616 | 0.097 |
| **RVPRS-FaST-LMM** | **0.414** | **0.629** | **0.097** |
| **RVPRS-HEELS** | **0.413** | **0.627** | **0.096** |
| OLS | 0.413 | 0.627 | 0.095 |
| IVW | 0.414 | 0.627 | 0.094 |
| SKAT-O | 0.414 | 0.626 | 0.095 |
| Burden | 0.546 | 0.340 | 0.028 |

RVPRS methods are competitive with individual-level baselines (Ridge, LASSO) despite using only summary statistics. Burden PRS underperforms due to its inability to estimate per-variant effects.

### Architecture comparison (quantitative, h2=0.15, matched LD)

Under **polygenic** architecture, RVPRS-FaST-LMM achieves RMSE 0.066 — matching OLS (0.065), the theoretical optimum when n >> p.

Under **sparse** architecture, LASSO/Elastic Net dominate (RMSE 0.187) as expected, but RVPRS-FaST-LMM (0.208) outperforms Ridge (0.221) due to better variance component estimation.

### LD mismatch stress test

| Condition | RVPRS-HEELS RMSE | RVPRS-FaST-LMM RMSE | Ridge RMSE |
|-----------|:-----------------:|:--------------------:|:----------:|
| Matched | 0.484 | 0.471 | 0.476 |
| Mismatched | 5.55 | 1.79 | 0.342 |

LD mismatch degrades summary-stat methods significantly. RVPRS-FaST-LMM is more robust than RVPRS-HEELS. Individual-level methods (Ridge) are completely unaffected. This validates the importance of reference panel quality for rare-variant summary-stat methods (consistent with COJO's finding that reference N >= 5,000 is recommended).

---

## Usage

### Full grid run (all 648 scenarios, 20 reps each)

```bash
PYTHONPATH=/data/home/seokhojeong/Meta_RVPRS \
python -m rvprs_pipeline.simulation_study \
  --outdir /path/to/sim_results \
  --n-reps 20 \
  --seed 42 \
  --verbose
```

### Quick test (single scenario)

```bash
PYTHONPATH=/data/home/seokhojeong/Meta_RVPRS \
python -m rvprs_pipeline.simulation_study \
  --outdir /path/to/test_out \
  --n-values 5000 \
  --p-values 40 \
  --architectures mixed \
  --h2-values 0.05 \
  --trait-types quantitative \
  --ld-mismatches matched \
  --cohort-values 1 \
  --n-reps 3
```

### Focused scenario subset

```bash
PYTHONPATH=/data/home/seokhojeong/Meta_RVPRS \
python -m rvprs_pipeline.simulation_study \
  --outdir /path/to/focused \
  --n-values 10000 50000 \
  --architectures sparse mixed \
  --h2-values 0.05 0.15 \
  --trait-types quantitative \
  --ld-mismatches matched \
  --cohort-values 1 \
  --n-reps 10 \
  --methods rvprs_heels rvprs_fastlmm burden ridge lasso
```

### Summary-only mode (skip oracle baselines)

```bash
--summary-only --methods rvprs_heels rvprs_fastlmm burden ivw skato_weights
```

---

## Files Modified / Created

| File | Action | Description |
|------|--------|-------------|
| `simulation_study.py` | **Created** | Comprehensive simulation study framework |
| `session_summary_20260416.md` | **Created** | This progress document |
| `simulation_study_demo.ipynb` | **Created** | Interactive notebook for running and visualizing results |

---

## Next Steps

1. **Run full grid** on a compute node (648 scenarios x 20 reps takes ~hours at N=200k)
2. **Add collapse analysis**: Evaluate methods after MAC-based variant collapsing
3. **Add RareEffect R bridge**: Enable `--include-rareeffect` for head-to-head with the R implementation
4. **Publication figures**: Power curves (method x architecture x N), LD mismatch degradation plots, h2 calibration QQ plots
