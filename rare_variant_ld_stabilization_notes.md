# Rare-Variant LD Stabilization: Theory, Existing Literature, and Practical Modeling Directions

## Reading guide

This note is written as a self-contained technical document for understanding why **linkage disequilibrium (LD)** becomes unstable for **rare variants**, what the existing literature has already established, and where new modeling ideas may fit.

The goal is not only to summarize papers, but to make the logic explicit enough that the reader can follow it line by line:

1. define LD carefully;
2. distinguish population LD from sample-estimated LD;
3. explain why rare variants create a fundamentally difficult estimation regime;
4. review major existing strategies for stabilization;
5. identify where graph/genealogy priors and higher-order moments might enter.

Throughout, I will use the phrase **LD stabilization** in a broad mathematical sense:

- reducing estimator bias;
- reducing estimator variance;
- improving conditioning of LD/covariance matrices;
- injecting external structure so that pairwise LD estimates are not driven entirely by sparse sample counts.

---

## 1. Problem setup: what does it mean to stabilize LD for rare variants?

### 1.1 Why do we care about LD matrices?

LD enters statistical genetics in several distinct ways.

1. **Single-variant GWAS fine-mapping and summary-statistics regression** use an LD matrix or correlation matrix to connect marginal association statistics to joint effects.
2. **Conditional analysis** uses LD to determine whether two association signals are independent.
3. **Polygenic modeling** often requires an estimate of local or genome-wide SNP correlation.
4. **Rare-variant aggregation tests** such as burden tests, SKAT-type tests, and meta-analysis frameworks often need score covariance matrices or related objects.
5. **Imputation and haplotype-based inference** implicitly or explicitly rely on multi-locus dependence structure.

Thus, LD is not merely a descriptive statistic. It is a core inferential object.

### 1.2 Why are rare variants different?

For common variants, empirical LD estimated from a sufficiently large ancestry-matched sample often behaves reasonably well. The cell counts in the corresponding contingency table are not too sparse, sample covariance is less noisy, and the pairwise correlation structure can often be estimated with acceptable precision.

For rare variants, the situation changes:

- minor allele counts are small;
- phase uncertainty increases;
- genotyping and imputation uncertainty become more consequential;
- ancestry mismatch matters more strongly;
- many observed pairwise correlations are dominated by sampling noise.

Therefore, when we speak about **rare-variant LD stabilization**, we mean: how can we infer a dependence structure that is more reliable than the raw empirical pairwise correlations computed from sparse, noisy data?

### 1.3 Three levels of stabilization

It is useful to separate stabilization into three levels.

#### Estimator-level stabilization
Correct the way pairwise LD is estimated.

Examples:

- bias correction under genotype uncertainty;
- shrinkage of noisy sample correlations;
- robust estimation under low minor allele count.

#### Matrix-level stabilization
Take a noisy matrix estimate and regularize it.

Examples:

- diagonal loading;
- off-diagonal shrinkage;
- banding, tapering, sparsification;
- conditioning or precision-matrix estimation.

#### Structure-level stabilization
Use biological or genealogical structure not visible in raw pairwise counts.

Examples:

- haplotype-copying models;
- genealogy-informed graphical models;
- region-level aggregation;
- higher-order dependence summaries.

This distinction will matter later, because much of the existing literature addresses the first two levels, whereas your ideas about graph priors and third moments live closer to the third level.

---

## 2. Basic definitions and notation

### 2.1 Two-locus setup

Consider two biallelic loci, indexed by 1 and 2.

Let the alleles at locus 1 be $A/a$ and at locus 2 be $B/b$.
Define haplotype frequencies

$$
p_{AB}, \quad p_{Ab}, \quad p_{aB}, \quad p_{ab},
$$

with

$$
p_{AB} + p_{Ab} + p_{aB} + p_{ab} = 1.
$$

Marginal allele frequencies are

$$
p_A = p_{AB} + p_{Ab}, \qquad p_B = p_{AB} + p_{aB}.
$$

### 2.2 The classical LD coefficient $D$

The classical two-locus LD coefficient is

$$
D = p_{AB} - p_A p_B.
$$

This is the difference between the observed frequency of the $AB$ haplotype and the frequency we would expect if the loci were independent.

If the two loci are in linkage equilibrium, then

$$
p_{AB} = p_A p_B,
$$

and therefore $D = 0$.

If $D \neq 0$, allelic states are statistically associated.

### 2.3 Standardized measures: $r$ and $r^2$

Because $D$ depends on allele frequencies, it is common to standardize it. Define

$$
r = \frac{D}{\sqrt{p_A(1-p_A)p_B(1-p_B)}},
$$

and

$$
r^2 = \frac{D^2}{p_A(1-p_A)p_B(1-p_B)}.
$$

Interpretations:

- $r$ is the correlation between allelic indicators at the two loci;
- $r^2$ is the squared correlation and is often used as a measure of tagging strength.

### 2.4 Genotype representation

For statistical genetics, it is often useful to represent a locus $j$ by an individual-level genotype variable $G_j$ taking values in $\{0,1,2\}$, where $G_j$ is the count of alternate alleles.

Define the centered genotype

$$
X_j = G_j - \mathbb{E}[G_j].
$$

Then the covariance between loci $j$ and $k$ is

$$
\Sigma_{jk} = \operatorname{Cov}(G_j, G_k) = \mathbb{E}[X_j X_k].
$$

The correlation matrix is

$$
R_{jk} = \frac{\Sigma_{jk}}{\sqrt{\Sigma_{jj}\Sigma_{kk}}}.
$$

This $R$ is often what GWAS software informally calls the LD matrix.

### 2.5 Population LD versus sample LD

This distinction is crucial.

- **Population LD** is the underlying dependence structure in the population.
- **Sample LD** is the estimate computed from finite data.

Given $n$ samples, the empirical covariance is

$$
\hat\Sigma = \frac{1}{n} X^\top X,
$$

where $X$ is the $n \times p$ matrix of centered genotypes.

The empirical correlation matrix is

$$
\hat R_{jk} = \frac{\hat\Sigma_{jk}}{\sqrt{\hat\Sigma_{jj}\hat\Sigma_{kk}}}.
$$

Even if the population LD is well defined, $\hat R$ may be extremely noisy when the alleles are rare.

---

## 3. Classical two-locus theory

### 3.1 Recombination reduces LD

In the simplest deterministic setting, if the recombination fraction between the two loci is $c$, then under random mating,

$$
D_{t+1} = (1-c) D_t.
$$

Thus,

$$
D_t = (1-c)^t D_0.
$$

This gives the standard intuition: recombination decays LD over generations.

### 3.2 Drift, mutation, and selection complicate the picture

In finite populations, random genetic drift generates stochastic fluctuations in haplotype frequencies. Mutation introduces new alleles, often initially at very low frequency. Selection changes the trajectories of allele frequencies and therefore changes the distribution of pairwise LD as well.

For common-variant intuition, one often thinks about the balance between recombination breaking LD and drift reintroducing it. But rare variants require more care because they have not had much time to diffuse through frequency space, and many observed pairs are in the very low-count regime.

### 3.3 Hill–Robertson and Hill–Weir perspective

Classical work by Hill, Robertson, and later Hill & Weir studied the behavior of LD statistics in finite populations. One of the most enduring insights is that pairwise LD statistics are themselves random variables with nontrivial variance, and that their expected behavior depends on population size, recombination, and sampling.

This matters because in rare-variant settings, the variance of the estimator is not a small nuisance term. It can dominate the observed value.

---

## 4. Why empirical LD breaks down for rare variants

### 4.1 Sparse contingency tables

Suppose the minor allele frequency at locus $j$ is $f_j$ and at locus $k$ is $f_k$.
Then the expected double-carrier frequency is roughly on the order of

$$
O(f_j f_k)
$$

when the loci are weakly associated.

If both are rare, say $f_j, f_k \ll 1$, then the expected number of double carriers in a sample of size $n$ is approximately

$$
n f_j f_k,
$$

which may be much smaller than 1.

That means the empirical covariance is often determined by whether we happened to observe 0, 1, or 2 double carriers.

### 4.2 Instability of the sample covariance

Write the empirical covariance as

$$
\hat\Sigma_{jk} = \frac{1}{n} \sum_{i=1}^n (G_{ij} - \bar G_j)(G_{ik} - \bar G_k).
$$

When $G_{ij}$ is almost always 0, the quantity is sparse and highly sensitive to a few carriers. The numerator is effectively driven by a very small subset of samples.

At the same time, the diagonal entries

$$
\hat\Sigma_{jj} \approx 2 f_j (1-f_j)
$$

become small for rare variants, so the correlation normalization can amplify noise further:

$$
\hat R_{jk} = \frac{\hat\Sigma_{jk}}{\sqrt{\hat\Sigma_{jj}\hat\Sigma_{kk}}}.
$$

Thus, even a small random fluctuation in $\hat\Sigma_{jk}$ may produce a deceptively large $\hat R_{jk}$.

### 4.3 “No observed LD” is not the same as “no LD”

A key conceptual mistake is to equate

- no co-observed rare alleles in the sample,

with

- genuine biological independence.

If the sample is not large enough to observe joint carriers, then the empirical estimate can be near zero even if the underlying haplotype or genealogy implies nontrivial dependence.

### 4.4 Phase uncertainty and calling uncertainty

For common variants, small phase errors often have modest effect. For rare variants, a single switch or calling error can substantially alter the estimated LD, because the relevant information is concentrated in very few chromosomes.

Similarly, if genotypes are uncertain and we use posterior means naively, covariance estimates can be attenuated or biased.

### 4.5 Ancestry mismatch

Rare variants are often more population-specific than common variants. Therefore, an external reference panel may poorly represent the haplotype background of the study sample. This means that reference-panel LD mismatch is often more severe for rare variants than for common variants.

---

## 5. Rare-variant LD theory in more detail

### 5.1 Why classical common-variant intuition is insufficient

Common-variant LD is often interpreted through equilibrium-style heuristics. Rare variants are different because they usually sit in an early phase of their stochastic life history.

A rare mutation is typically young. Its distribution across haplotypes reflects recent lineage history, not a long-run equilibrium. Therefore, the observed LD between two rare variants depends strongly on their ages, shared genealogy, and whether they arose on related lineages.

### 5.2 Frequency-conditioned LD

A key insight from recent rare-variant LD theory is that LD should be considered conditional on present-day frequencies. In other words, the distribution of $D$, $r$, or related quantities is very different when both alleles are forced to be rare.

This conditioning changes the effective geometry of the problem. You are no longer averaging over all pairs of loci. You are averaging over a selected subset of pairs whose allele counts are tiny.

### 5.3 Rare mutation limit

In the rare mutation limit, one can often simplify the dynamics by tracking the branching behavior of the lineages carrying the mutations. This yields intuition of the following sort:

- if two rare variants arose on independent lineages, co-occurrence is extremely unlikely;
- if one arose on a background carrying the other, strong positive association may be observed;
- if selection acts differently on the two mutations, the distribution of LD changes in a frequency-dependent way.

This is one reason rare-variant LD is not just “common-variant LD but noisier.” The underlying theoretical regime is genuinely different.

### 5.4 Practical implication

A stabilization method that ignores frequency conditioning may systematically mischaracterize rare-variant LD. A method that treats all pairwise entries identically may be mathematically convenient but biologically unrealistic.

This motivates adaptive procedures in which shrinkage strength, prior variance, or graph connectivity depends on allele count, local genealogy, or local recombination context.

---

## 6. From population covariance to sample LD matrices

### 6.1 The high-dimensional setting

In modern genetics, the number of variants $p$ in a locus or region may be comparable to or larger than the sample size $n$. Then the sample covariance matrix

$$
\hat\Sigma = \frac{1}{n}X^\top X
$$

may be ill-conditioned or singular.

Even when $p < n$, if many variants are rare, effective information is much smaller than $n$, because most rows contribute no information to many off-diagonal pairs.

### 6.2 What goes wrong numerically?

A raw empirical LD matrix may suffer from:

- small or zero eigenvalues;
- unstable inverse or pseudo-inverse;
- noisy off-diagonal entries;
- poor portability across cohorts.

When such a matrix is used in regression, fine-mapping, or conditional analysis, the downstream estimates can become unstable.

### 6.3 Misspecification in summary-statistics models

Suppose marginal GWAS statistics $z$ are modeled with an LD matrix $R$ under some approximate likelihood.
If the matrix supplied to the model is actually $\tilde R \neq R$, then effect estimates, posterior inclusion probabilities, credible sets, and conditional analyses can all be distorted.

For rare variants, this misspecification may come less from classical estimation noise and more from systematic under-observation of pairwise dependence.

---

## 7. What LD stabilization can mean mathematically

Here are several mathematically distinct ideas that all count as stabilization.

### 7.1 Bias correction

If the estimator is biased because of genotype uncertainty or posterior mean compression, one can replace the naive estimator by a corrected one.

Symbolically,

$$
\hat R^{\text{corr}} = \mathcal{C}(\hat R^{\text{naive}}).
$$

### 7.2 Shrinkage

A general shrinkage estimator takes the form

$$
\hat\Sigma^{\text{shrink}} = (1-\lambda)\hat\Sigma + \lambda T,
$$

where $T$ is a target matrix.

Common choices:

- $T = \operatorname{diag}(\hat\Sigma)$;
- $T = I$ after standardization;
- local banded or recombination-informed targets.

The effect is to reduce variance and improve conditioning.

### 7.3 Banding or tapering

If distant pairs are expected to have small LD, one can damp them as a function of genomic distance:

$$
\hat\Sigma^{\text{taper}}_{jk} = w_{jk} \hat\Sigma_{jk},
$$

where $w_{jk}$ decreases with distance.

### 7.4 Sparse precision modeling

Instead of modeling covariance directly, one can model the precision matrix

$$
\Theta = \Sigma^{-1}.
$$

Sparsity in $\Theta$ corresponds to conditional independence structure. In statistical genetics, sparse precision ideas become appealing when dense LD matrices are too expensive or too noisy.

### 7.5 Structure-aware priors

Suppose we have a graph $G$ or some genealogy-derived object. We might define a prior mean matrix $M(G)$ and shrink the empirical LD toward it:

$$
\hat\Sigma^{\text{post}} = (1-\lambda)\hat\Sigma + \lambda M(G).
$$

This is closer to the direction you are considering.

### 7.6 Region-level aggregation

Instead of stabilizing pairwise SNP-by-SNP LD entries, define region-level quantities that naturally average over rare variants and therefore have lower variance.

This leads to statistics such as cumulative LD (cLD).

---

## 8. Bias-corrected LD estimation under genotype uncertainty

### 8.1 The problem

Suppose the genotype at locus $j$ for individual $i$ is not observed exactly, but represented by a posterior distribution or genotype likelihood. A common shortcut is to replace the latent genotype by its posterior mean:

$$
\tilde G_{ij} = \mathbb{E}[G_{ij} \mid \text{data}].
$$

One then computes covariance or correlation using $\tilde G_{ij}$ as if it were the true genotype.

But this generally induces bias, because

$$
\operatorname{Cov}(\mathbb{E}[G_j \mid D], \mathbb{E}[G_k \mid D])
\neq \operatorname{Cov}(G_j, G_k).
$$

The posterior mean smooths over uncertainty and tends to attenuate variation.

### 8.2 Moment-based correction intuition

Using the law of total covariance,

$$
\operatorname{Cov}(G_j, G_k)
= \operatorname{Cov}(\mathbb{E}[G_j \mid D], \mathbb{E}[G_k \mid D]) + \mathbb{E}[\operatorname{Cov}(G_j, G_k \mid D)].
$$

If the conditional covariance term is ignored, the estimate is biased. The bias-corrected approach reintroduces the missing uncertainty contribution using tractable moment adjustments.

### 8.3 Why this matters especially for rare variants

Rare variants are exactly where genotype uncertainty hurts most:

- low read depth;
- uncertain imputation dosage;
- sparse carrier counts;
- phase uncertainty.

Therefore, a practical stabilization pipeline should often start with an uncertainty-aware estimator before any matrix regularization is attempted.

### 8.4 Conceptual takeaway

Not all apparent instability is biological. Some is pure measurement/estimation error. It makes little sense to design sophisticated graph priors if the base pairwise estimator is already biased because of genotype uncertainty.

---

## 9. Shrinkage and regularization of LD matrices

### 9.1 Generic shrinkage

A shrinkage covariance estimator has the form

$$
\hat\Sigma_\lambda = (1-\lambda)\hat\Sigma + \lambda T,
$$

with $0 \leq \lambda \leq 1$.

As $\lambda$ increases, the estimate becomes less variable but potentially more biased. The art is to choose a target $T$ and strength $\lambda$ that reduce noise while preserving real structure.

### 9.2 Off-diagonal shrinkage

One particularly natural choice for LD is to shrink off-diagonal entries toward zero while preserving marginal variances:

$$
\hat\Sigma^{\text{off-shrink}}_{jk} =
\begin{cases}
\hat\Sigma_{jj}, & j = k, \\
(1-\lambda_{jk})\hat\Sigma_{jk}, & j \neq k.
\end{cases}
$$

This is reasonable when many observed off-diagonal terms are dominated by noise.

### 9.3 Why shrinkage helps

It improves:

- conditioning of the matrix;
- invertibility or approximate invertibility;
- robustness in downstream regression;
- portability across samples.

### 9.4 What is missing for rare variants?

Standard shrinkage methods typically do not explicitly model the fact that a pair with minor allele counts $(2,2)$ should be treated differently from a pair with counts $(500,700)$.

This suggests rare-variant-aware shrinkage strengths of the form

$$
\lambda_{jk} = \lambda(\text{MAC}_j, \text{MAC}_k, d_{jk}, \text{uncertainty}, \text{ancestry match}),
$$

where $d_{jk}$ is genomic distance.

### 9.5 A possible extension

A principled rare-variant extension could combine:

1. uncertainty-aware pairwise estimation;
2. allele-count-adaptive shrinkage;
3. graph-informed prior mean.

Then the final estimate would look like

$$
\hat\Sigma^{\star}_{jk} = (1-\lambda_{jk})\hat\Sigma^{\text{corr}}_{jk} + \lambda_{jk} M_{jk}(G),
$$

where $M(G)$ is derived from some graph or genealogy object.

---

## 10. Sparse graph and precision-based representations

### 10.1 Why switch from covariance to precision?

A covariance matrix tells us pairwise marginal dependence. A precision matrix

$$
\Theta = \Sigma^{-1}
$$

tells us conditional dependence. In many settings, conditional dependence is sparser and more stable to model.

If a dense LD matrix is hard to estimate reliably, it may be easier to infer a sparse graphical representation that captures the dominant conditional structure.

### 10.2 Gaussian graphical model analogy

In a Gaussian graphical model, an edge between nodes $j$ and $k$ exists when

$$
\Theta_{jk} \neq 0.
$$

Although genotype data are not Gaussian, this analogy is still useful. A sparse precision object may serve as an efficient surrogate for a dense correlation matrix.

### 10.3 LD graphical models (LDGM)

Recent work introduced **LD graphical models**, which construct sparse graph-based representations of LD derived from genome-wide genealogies. The main idea is not that the true LD matrix is sparse, but that a sparse graph may approximate its algebraic behavior very efficiently.

This is important for your question because it shows that:

- graph-based LD modeling is already a serious direction in the literature;
- genealogy can be used as a structural scaffold for LD;
- sparse representations can outperform naive dense empirical matrices computationally and sometimes statistically.

### 10.4 Limitation for rare variants

Existing LDGM work is more developed for common variants than for the ultra-rare regime. So the conceptual template exists, but the rare-variant extension remains open territory.

---

## 11. Genealogy-based modeling: ARGs, local trees, and tree sequences

### 11.1 ARG intuition

An **ancestral recombination graph (ARG)** describes how sampled chromosomes coalesce backward in time while allowing recombination. Because recombination changes local ancestry along the genome, different genomic intervals may have different local trees.

### 11.2 Why genealogy matters for LD

LD is fundamentally a statement about shared ancestry and recombination. If two loci tend to share genealogy, they may exhibit stronger dependence. If they often separate by recombination, LD decays.

This is why genealogical interpretations of LD are so powerful: they explain LD as a consequence of correlated coalescent histories.

### 11.3 Local trees and tree sequences

In practice, one often works with a sequence of local trees along the genome rather than an explicit full ARG object. A **tree sequence** is an efficient representation of these local genealogies.

### 11.4 Why genealogy is attractive for rare variants

Rare variants are often recent and lineage-specific. Their dependence may therefore be more naturally expressed relative to the local genealogy than relative to a raw pairwise sample correlation table.

This suggests a key conceptual move:

> Instead of estimating rare-variant LD purely from observed co-carriers, estimate it by combining sparse observed data with a prior derived from local genealogy.

### 11.5 Link to Li–Stephens-style copying models

The Li–Stephens model explains haplotypes through copying from existing haplotypes with occasional recombination. This is not a full ARG, but it captures the same broad intuition: local dependence is mediated by mosaic copying along ancestral backgrounds.

For stabilization, this means one can borrow information from the latent haplotype structure even when a rare pair itself has almost no observed co-occurrences.

---

## 12. Linkage or genealogy graph as a prior for rare-variant LD

This section is closest to the first of your proposed directions.

### 12.1 Motivation

Suppose the empirical pairwise LD entry $\hat\Sigma_{jk}$ is very noisy because the rare alleles appear in too few individuals. But suppose we have a graph $G$ on variants or haplotype states, representing local genealogical or copying relationships.

Then it may be sensible to define a prior mean $M(G)$ and shrink the empirical estimate toward it.

### 12.2 Basic shrinkage-with-prior form

A generic formulation is

$$
\hat\Sigma^{\text{stab}} = (1-\Lambda) \odot \hat\Sigma^{\text{corr}} + \Lambda \odot M(G),
$$

where:

- $\odot$ denotes elementwise multiplication;
- $\Lambda = (\lambda_{jk})$ is an entry-specific shrinkage matrix;
- $\hat\Sigma^{\text{corr}}$ is a corrected sample estimator;
- $M(G)$ is a graph-derived prior mean.

Entrywise,

$$
\hat\Sigma^{\text{stab}}_{jk} = (1-\lambda_{jk})\hat\Sigma^{\text{corr}}_{jk} + \lambda_{jk} M_{jk}(G).
$$

### 12.3 How could $M(G)$ be defined?

Several possibilities exist.

#### Path-based prior

If $d_G(j,k)$ is a graph distance, define

$$
M_{jk}(G) = \alpha \exp(-\beta d_G(j,k)).
$$

This says that nearby nodes in the genealogy or linkage graph should have higher prior covariance.

#### Laplacian smoothness prior

Let $L_G$ be the graph Laplacian. One may penalize rough covariance structure across the graph via an optimization such as

$$
\min_\Sigma \; \|\Sigma - \hat\Sigma^{\text{corr}}\|_F^2 + \tau \, \operatorname{tr}(\Sigma L_G \Sigma^\top),
$$

subject to positive semidefiniteness.

The second term encourages covariance structure that is smooth with respect to the graph.

#### Precision prior

Instead of placing a prior on $\Sigma$, one could place a sparsity or smoothness prior on

$$
\Theta = \Sigma^{-1},
$$

and use graph adjacency to determine which off-diagonal precision terms are allowed.

### 12.4 Why this might help rare variants

For ultra-rare pairs, the raw data contain almost no information about $\Sigma_{jk}$. A graph prior lets the estimate borrow information from neighboring variants or shared haplotype backgrounds.

This is conceptually similar to empirical Bayes smoothing: low-information entries get shrunk more strongly toward a structurally informed target.

### 12.5 Caveat

A graph prior can help only if the graph is itself informative and ancestry-matched. A mismatched graph may inject systematic bias.

---

## 13. Beyond pairwise LD: third moments and higher-order dependence

This is closest to your second proposed direction.

### 13.1 Why pairwise correlation may be insufficient

Pairwise covariance captures only second-order dependence. But haplotype structure is fundamentally multi-locus. In principle, two loci can appear weakly correlated marginally, while higher-order dependence involving a third locus reveals shared structure.

### 13.2 Third central moments

For centered genotype variables $X_j$, define the third central moment

$$
\mu_{jkl}^{(3)} = \mathbb{E}[X_j X_k X_l].
$$

This measures triadic dependence. If all variables were jointly Gaussian with mean zero, this quantity would vanish. But genotype data are not Gaussian, and tri-locus interactions or haplotype constraints may generate nonzero values.

### 13.3 Cumulant perspective

The third-order cumulant is the same as the third central moment for centered variables. Higher-order cumulants measure structure not explained by lower-order moments.

In that sense, pairwise LD explains only a slice of multi-locus dependence.

### 13.4 Why direct estimation is hard for rare variants

If rare-rare pairwise counts are sparse, then rare-rare-rare triple counts are even sparser. The naive plug-in estimate of $\mu_{jkl}^{(3)}$ may have enormous variance.

In a sample of size $n$, the number of individuals carrying all three rare alleles may be near zero for the overwhelming majority of triplets.

Thus, direct third-moment estimation is statistically fragile.

### 13.5 More realistic role for third moments

Because direct estimation is so hard, the best use of third moments may be indirect.

Possible roles:

1. **Hyperparameter learning**: use aggregated third-moment summaries to decide where pairwise shrinkage should be weaker or stronger.
2. **Motif detection**: identify local haplotype motifs or graph substructures where pairwise LD is under-estimated.
3. **Region-level diagnostics**: determine whether a region exhibits multi-locus dependence beyond what pairwise LD explains.

### 13.6 A conceptual model

One could imagine

$$
\lambda_{jk} = f\left(\text{MAC}_j, \text{MAC}_k, d_{jk}, U_{jk}\right),
$$

where $U_{jk}$ is a summary derived from third-order statistics involving loci near $j$ and $k$.

For example,

$$
U_{jk} = \sum_{l \in \mathcal{N}(j,k)} w_l \, |\hat\mu^{(3)}_{jkl}|.
$$

Then higher-order information does not replace pairwise LD; it modulates the stabilization procedure.

This is much more plausible than trying to estimate and use a full third-order tensor directly.

---

## 14. Region-level alternatives: cumulative LD and other aggregation ideas

### 14.1 Why move to region-level statistics?

If single-entry pairwise LD estimates are intrinsically unstable, perhaps the object being estimated should be changed.

This is the logic behind region-level rare-variant LD statistics.

### 14.2 Cumulative LD (cLD)

Recent work proposed **cumulative LD (cLD)** as a statistic that captures LD between genomic regions containing rare variants. The goal is to obtain a dependence measure whose variance is lower and whose biological interpretation remains meaningful.

The core intuition is simple:

- SNP-level rare-rare LD is too sparse;
- aggregating across multiple variants increases effective information;
- the resulting statistic is more stable than pairwise $r$ or $r^2$.

### 14.3 General region-level formulation

Suppose $A$ and $B$ are two sets of variants. Define collapsed genotype burdens

$$
S_A = \sum_{j \in A} w_j X_j, \qquad S_B = \sum_{k \in B} v_k X_k.
$$

Then a region-level dependence measure is naturally

$$
\operatorname{Cov}(S_A, S_B) = \sum_{j \in A} \sum_{k \in B} w_j v_k \Sigma_{jk}.
$$

Even if each individual $\Sigma_{jk}$ is noisy, the aggregate may be much more stable.

### 14.4 Why this matters for your problem

If your actual downstream goal is gene-level testing, pathway testing, or regional interpretation, then stabilizing a region-level covariance operator may be more useful than stabilizing every matrix entry.

This is an important design question:

> Do we truly need a dense pairwise LD matrix for rare variants, or do we need a stable dependence object for downstream inference?

These are not the same question.

---

## 15. Connection to rare-variant association testing

### 15.1 Score statistics and covariance

For many rare-variant tests, the starting point is a vector of single-variant score statistics $U$ with covariance matrix $V$.

A quadratic-form test often takes the form

$$
Q = U^\top W U,
$$

where the null distribution depends on the covariance structure of $U$.

Thus, even when the analysis is not phrased directly in terms of an LD matrix, a covariance object remains central.

### 15.2 Burden and SKAT-style tests

Burden tests effectively aggregate variants before testing. SKAT-type methods preserve more local heterogeneity but still require covariance information.

In summary-statistics meta-analysis settings, efficient storage and use of these covariance matrices become a major practical problem.

### 15.3 Implication

The rare-variant testing literature reinforces a broad lesson:

- dense exact pairwise LD may be too expensive or unstable;
- sparse, structured, or aggregated covariance objects are often sufficient and preferable.

This aligns well with graph-based stabilization and region-level alternatives.

---

## 16. Practical sources of error in real data

### 16.1 Small sample size

The most obvious problem. Rare-variant double carriers are too few.

### 16.2 Low minor allele count

A variant can have acceptable nominal MAF in principle yet still have low effective information in the analyzed subgroup.

### 16.3 Genotype calling uncertainty

Posterior means may attenuate variance and covariance.

### 16.4 Phasing error

Rare-variant LD is strongly haplotype dependent, so phase errors matter more.

### 16.5 Imputation error

Reference mismatch and low imputation quality can destroy fine-scale rare-variant dependence.

### 16.6 Ancestry mismatch

Rare variants are often ancestry-specific and localized to particular genealogical backgrounds.

### 16.7 Local recombination heterogeneity

Hotspots and local structure can make distance-based smoothing too crude.

### 16.8 Structural variation and complex loci

Inversion regions, segmental duplications, and structural polymorphisms can induce unusual patterns not captured by simple pairwise decay assumptions.

---

## 17. A unified view of stabilization

It is useful to summarize the whole landscape in one table-like conceptual framework.

### 17.1 Level 1: estimator correction

Goal: improve raw pairwise estimates.

Tools:

- genotype-uncertainty correction;
- better phase-aware estimation;
- adaptive variance estimates.

### 17.2 Level 2: matrix regularization

Goal: improve conditioning and reduce noise.

Tools:

- shrinkage;
- tapering;
- sparsification;
- precision-based methods.

### 17.3 Level 3: structural priors or alternative objects

Goal: encode biology not visible in sparse sample co-counts.

Tools:

- genealogy graph prior;
- LDGM-like sparse representations;
- higher-order summaries;
- region-level dependence statistics.

### 17.4 Most plausible practical combination

For rare variants, the most sensible workflow is likely not a single magical estimator, but a stack:

1. **correct the raw estimator**;
2. **regularize the matrix**;
3. **borrow structure from graph/genealogy or aggregate regionally**.

---

## 18. Candidate modeling frameworks for the current use case

Here are four concrete ways to proceed.

### Framework A: uncertainty-aware sample LD plus adaptive shrinkage

1. Estimate pairwise covariance using genotype-uncertainty correction.
2. Apply entrywise shrinkage with strength depending on minor allele counts and distance.
3. Use the result as the stabilized LD matrix.

A stylized estimator:

$$
\hat\Sigma^{A}_{jk} = (1-\lambda_{jk})\hat\Sigma^{\text{corr}}_{jk}, \qquad j \neq k.
$$

Pros:

- simple;
- easy baseline;
- likely publishable as a careful empirical study if done well.

Cons:

- limited biological structure;
- still pairwise only.

### Framework B: genealogy-informed graph prior

1. Build or import a local linkage/genealogy graph.
2. Define a prior mean or smoothness operator from the graph.
3. Combine corrected sample LD with graph prior.

Estimator:

$$
\hat\Sigma^{B}_{jk} = (1-\lambda_{jk})\hat\Sigma^{\text{corr}}_{jk} + \lambda_{jk}M_{jk}(G).
$$

Pros:

- biologically interpretable;
- directly aligned with recent LDGM/ARG ideas;
- especially attractive when sample information is weak.

Cons:

- depends on graph quality and ancestry matching;
- needs careful validation.

### Framework C: pairwise LD with higher-order-informed hyperparameters

1. Compute pairwise corrected LD.
2. Estimate local summary measures of tri-locus dependence.
3. Use those summaries to modulate shrinkage or graph weights.

This avoids direct use of a full third-order tensor while still incorporating information beyond pairwise moments.

Pros:

- scientifically novel;
- potentially captures local haplotype motifs not visible in pairwise statistics.

Cons:

- statistically delicate;
- requires careful design to avoid noise amplification.

### Framework D: region-level operator instead of SNP-by-SNP LD

1. Accept that rare-variant SNP-level LD is intrinsically unstable.
2. Define region-level covariance objects or cumulative LD.
3. Use those for downstream testing or interpretation.

Pros:

- stable;
- natural for gene-based tests;
- close to existing rare-variant practice.

Cons:

- loses some pairwise resolution;
- less directly reusable in single-variant fine-mapping style workflows.

---

## 19. What the literature has already answered, and what remains open

### 19.1 What is well established

The literature already strongly supports the following claims.

1. **LD can be interpreted genealogically**, not only as an empirical correlation table.
2. **Rare-variant LD lives in a different theoretical regime** from common-variant LD.
3. **Raw empirical LD is unstable for rare variants**, especially under low counts or genotype uncertainty.
4. **Bias-corrected estimation under uncertainty is important**.
5. **Shrinkage and sparse representations help** for LD matrices more generally.
6. **Region-level alternatives can be much more stable** for rare variants.

### 19.2 What appears still open

These questions remain comparatively open and promising.

1. A principled **rare-variant-specific pairwise LD shrinkage estimator**.
2. A **genealogy-informed prior** specifically designed for rare-variant LD stabilization.
3. A rigorous way to use **third moments or tri-locus structure** without incurring prohibitive variance.
4. A comparison of **pairwise stabilization versus region-level replacement** for downstream tasks.

These open points are exactly why your directions look worthwhile.

---

## 20. Recommended reading order

If you want the fastest route to a strong conceptual foundation, I would read in this order.

1. **Li & Stephens (2003)** for multi-locus LD modeling through haplotype copying.
2. **McVean (2002)** for genealogical interpretation of LD.
3. **Good (2022)** for rare-mutation LD theory.
4. **Gerard et al. (2021)** for bias-corrected LD estimation under genotype uncertainty.
5. **Zhu & Stephens (2017)** for practical shrinkage/regularization context in summary-statistics analysis.
6. **O'Connor et al. (2023)** for LD graphical models and genealogy-informed sparse structure.
7. **Wang et al. (2023)** for cumulative LD as a rare-variant-stable alternative.
8. **MetaSTAAR (2022)** and related rare-variant meta-analysis papers to see how covariance objects are used in practice.

---

## 21. Concluding perspective

The central message of this note is the following.

Rare-variant LD instability is not just a numerical annoyance. It is the visible symptom of a deeper mismatch between:

- the object we want, namely a biologically meaningful multi-locus dependence structure;
- and the object we usually estimate, namely a raw pairwise sample correlation matrix.

When alleles are rare, the raw matrix is often under-informed. Stabilization therefore has to do more than denoise. It has to **replace missing information with principled structure**.

There are three increasingly ambitious ways to do that.

1. **Correct the raw estimator** so that technical uncertainty is not mistaken for biology.
2. **Regularize the matrix** so that downstream inference is numerically stable.
3. **Inject biological structure**, such as genealogy, graph organization, or region-level aggregation, so that the estimator reflects the latent haplotype architecture even when direct co-count evidence is weak.

Of your two original ideas, the current literature suggests:

- **graph or genealogy priors** are the more immediately grounded and likely tractable route;
- **third moments** are intellectually interesting, but probably most useful as an auxiliary signal rather than as a direct replacement for pairwise LD.

That does not mean the third-moment direction is weak. It means the best first implementation may be:

> use higher-order dependence to learn where and how much to shrink, rather than trying to estimate a full higher-order dependence object directly.

In short, the likely progression is:

$$
\text{corrected pairwise LD} \; + \; \text{adaptive shrinkage} \; + \; \text{genealogy-informed prior} 
\quad \Longrightarrow \quad \text{practical rare-variant LD stabilization}.
$$

---

## References

Below is a focused reference list for the topics discussed above.

### Foundational LD and genealogy

1. **Li, N. & Stephens, M. (2003).** Modeling Linkage Disequilibrium and Identifying Recombination Hotspots Using Single-Nucleotide Polymorphism Data. *Genetics*, 165(4), 2213-2233.  
   DOI: 10.1093/genetics/165.4.2213  
   Link: <https://academic.oup.com/genetics/article-abstract/165/4/2213/6050566>

2. **McVean, G. A. T. (2002).** A Genealogical Interpretation of Linkage Disequilibrium. *Genetics*, 162(2), 987-991.  
   DOI: 10.1093/genetics/162.2.987  
   Link: <https://academic.oup.com/genetics/article/162/2/987/6049967>

3. **Hill, W. G. & Robertson, A. (1968).** Linkage disequilibrium in finite populations. *Theoretical and Applied Genetics*, 38, 226-231.  
   Link: <https://link.springer.com/content/pdf/10.1007/BF01245622.pdf>

### Rare-variant LD theory

4. **Good, B. H. (2022).** Linkage disequilibrium between rare mutations. *Genetics*, 220(4), iyac004.  
   DOI: 10.1093/genetics/iyac004  
   Link: <https://academic.oup.com/genetics/article/220/4/iyac004/6503502>

### Bias-corrected LD estimation

5. **Gerard, D. et al. (2021).** Scalable bias-corrected linkage disequilibrium estimation under genotype uncertainty. *Heredity*, 127, 357-372.  
   Link: <https://www.nature.com/articles/s41437-021-00462-5>

### Shrinkage and summary-statistics context

6. **Zhu, X. & Stephens, M. (2017).** Bayesian large-scale multiple regression with summary statistics from genome-wide association studies. *Annals of Applied Statistics*, 11(3), 1561-1592.  
   Link: <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-3/Bayesian-large-scale-multiple-regression-with-summary-statistics-from-genome/10.1214/17-AOAS1046.full>

### Sparse graph / genealogy-informed LD modeling

7. **O'Connor, L. J. et al. (2023).** Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies. *Nature Genetics*, 55, 1494-1502.  
   Link: <https://www.nature.com/articles/s41588-023-01487-8>

8. **LDGM data resource (2023).** Data from Extremely sparse models of linkage disequilibrium in ancestrally diverse association studies. Zenodo.  
   Link: <https://zenodo.org/records/8157131>

### Region-level rare-variant LD

9. **Wang, X. et al. (2023).** cLD: Rare-variant linkage disequilibrium between genomic regions identifies novel molecular interactions. *PLOS Genetics*, 19(12), e1011074.  
   Link: <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011074>

### Rare-variant summary covariance and meta-analysis

10. **Li, X. et al. (2022).** Powerful, scalable and resource-efficient meta-analysis of rare variant associations in large whole-genome sequencing studies. *Nature Genetics*, 55, 154-164. [MetaSTAAR]  
    Link: <https://www.nature.com/articles/s41588-022-01225-6>

---

## Optional next expansions

If this note is extended further, the next sections I would add are:

1. a worked toy example showing how low minor allele count inflates variance of $\hat R_{jk}$;
2. an explicit derivation of uncertainty-induced attenuation bias;
3. a more formal Bayesian model for graph-prior shrinkage;
4. a simulation design section comparing raw LD, shrinkage LD, graph-prior LD, and cLD-based alternatives.

