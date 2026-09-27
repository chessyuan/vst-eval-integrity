# VST-Eval-Integrity: Competitor-Aware Evaluation Integrity Benchmark for Virtual Spatial Transcriptomics

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Article Track](https://img.shields.io/badge/Nature%20Methods-Analysis-red.svg)](https://www.nature.com/nmeth/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**VST-Eval-Integrity** is a fully reproducible, modular computational benchmark suite designed to audit data leakage, coordinate shortcuts, platform resolution divides, and foundation model scaling in **Virtual Spatial Transcriptomics (VST)**. 

Targeted for submission to ***Nature Methods*** as an **Analysis** article (in accordance with journal guidelines for comprehensive multi-method benchmarks), this repository operationalizes the author-proposed **VST-Eval Checklist** across 76 clinical tissue slides, 46 unique biological donors, 5 benchmarked foundation and classical architectures, and 4,632 empirical spatial permutation tests.

---

## Repository Structure

```
VST-Eval-Integrity/
├── configs/
│   └── benchmark_config.yaml         # Central execution parameters, thresholds, and cohort paths
├── donor_audit/
│   ├── audit_donor_provenance.py     # Biological donor verification and serial section quarantine
│   └── DONOR_PROVENANCE_MASTER.tsv   # Harmonized patient metadata across 76 slides and 46 donors
├── splits/
│   ├── leakage_ladder.py             # 5-level leakage ladder audit and patient-paired EOR engine
│   ├── leakage_ladder_metrics.tsv    # Evaluation performance from Level 0 to Level 6
│   └── patient_paired_eor.tsv        # Paired evaluation optimism metrics across 24 CCRCC patients
├── coordinate_nulls/
│   ├── coordinate_models.py          # Spatial coordinate baselines (C0-C5) and 2x2 factorial audit
│   ├── coordinate_null_suite.tsv     # Per-patient performance across 6 coordinate representations
│   └── factorial_image_coord.tsv     # Factorial interaction table testing feature fusion
├── spatial_nulls/
│   ├── spatial_permutation.py        # Variogram-preserving spatial permutation null engine
│   └── spatial_autocorrelation_null.tsv # Slide-level morphology retention and Moran's I metrics
├── equivalence_tests/
│   ├── tost_engine.py                # Two One-Sided Tests (TOST) and mixed-effects scaling engine
│   ├── tost_equivalence_tests.tsv   # Formal equivalence statistics at +/- 0.003 & +/- 0.005
│   ├── model_scaling_mixed_effects.tsv # Linear mixed-effects model regression table
│   └── pathway_vs_gene_harmonized_paired_test.tsv # Harmonized pathway fidelity paired tests
├── support_harmonization/
│   ├── hexagonal_binning.py          # Point-spread convolution ablation & exploratory ANOVA
│   ├── support_harmonization_ablation.tsv # 4-arm cross-platform resolution ablation
│   └── variance_decomposition_platform.tsv # ANOVA variance partitioning with 95% bootstrap CIs
├── uncertainty/
│   ├── selective_prediction.py       # Epistemic uncertainty & selective prediction retention
│   └── table5_uncertainty_failure_detection.tsv # AUROC and selective prediction metrics
├── model_adapters/
│   ├── adapters.py                   # Unified inference adapters (DeepSpot, HisToGene, ST-Net, etc.)
│   └── MODEL_STATUS_20260926.tsv     # Comprehensive audit of 13 published model repositories
├── reproduce_figures/
│   ├── generate_benchmark_figures.py # Script generating 6 publication figures (300 DPI)
│   ├── fig1_leakage_ladder_eor.png   # Figure 1: Leakage ladder & patient-paired EOR
│   ├── fig2_coordinate_null_tost.png # Figure 2: Coordinate nulls & TOST equivalence
│   ├── fig3_spatial_permutation_nulls.png # Figure 3: Spatial nulls & morphology retention
│   ├── fig4_support_harmonization.png# Figure 4: Support convolution & ANOVA decomposition
│   ├── fig5_model_scaling_deficit.png# Figure 5: Model scaling & structural confounding
│   └── fig6_pathway_harmonized_gain.png # Figure 6: Pathway fidelity & uncertainty retention
├── environment.yml                   # Complete Conda environment specification
├── Snakefile                         # Master Snakemake workflow definition
└── README.md                         # Benchmark documentation and replication guide
```

---

## Installation & Quick Start

### 1. Conda Environment Setup
Create and activate the standardized Conda environment:
```bash
conda env create -f environment.yml
conda activate vst-eval
```

### 2. End-to-End Execution via Snakemake
To run all verification stages, re-compute audit statistics, and regenerate all publication figures:
```bash
snakemake benchmark --cores 8
```

### 3. Standalone Script Execution
Each module can be executed independently from the command line:
```bash
# 1. Biological Donor Provenance Audit
python donor_audit/audit_donor_provenance.py

# 2. Evaluation Leakage Ladder & Patient-Paired EOR
python splits/leakage_ladder.py

# 3. Coordinate Null Baselines & Factorial Analysis
python coordinate_nulls/coordinate_models.py

# 4. Variogram-Preserving Spatial Permutation Nulls
python spatial_nulls/spatial_permutation.py

# 5. TOST Equivalence Testing & Parameter Scaling Mixed-Effects
python equivalence_tests/tost_engine.py

# 6. Support Harmonization & Exploratory Variance Partitioning
python support_harmonization/hexagonal_binning.py

# 7. Epistemic Uncertainty & Selective Prediction
python uncertainty/selective_prediction.py

# 8. Unified Model Adapters Test
python model_adapters/adapters.py

# 9. Publication Figure Reproduction
python reproduce_figures/generate_benchmark_figures.py
```

---

## Core Scientific Findings & Statistical Calibrations

### 1. Spatial Autocorrelation & Morphology Retention (Denominator-Grounded)
- **Cohort Scope**: Evaluated across 24 CCRCC patient slides $\times$ 193 target genes = **4,632 total gene–slide combinations**.
- **Empirical Significance**: Across 1,000 variogram-preserving spatial permutations per slide, exactly **373 of 4,632 tests (8.05%)** exhibited statistically significant morphology-associated predictive signal exceeding both spatial nulls (empirical FDR $< 0.05$) and coordinate baselines.
- **Slide Distribution**: Reported strictly as:
  $$\text{mean } 15.54 \text{ of } 193 \text{ genes per slide } (8.05\% \text{ of } 4,632 \text{ tested combinations}; 95\% \text{ CI: } [7.65, 23.44], \text{SD} = 19.73)$$

### 2. Negative Association in Foundation Model Scaling (Causal Purge)
- **Mixed-Effects Model**: Across 204 patient-level model evaluations, parameter scaling followed:
  $$\text{Gene PCC} = 0.0777 - 0.0430 \times \log_{10}(\text{Parameters}) + \text{Generalization Level} + (1 | \text{Patient})$$
  with $\beta = -0.0430$ (SE = 0.0034, 95% CI: [-0.0496, -0.0364], $p = 1.37 \times 10^{-37}$).
- **Scientific Calibration**: Causal claims (*"proves scaling hurts generalization"*) are rejected. Parameter scale co-varies with architecture family (CNN vs. ViT), pretraining data distribution, and training objectives (**structural confounding**).

### 3. Source-Domain / Acquisition-Domain Predictability
- **Confounding Breakdown**: High multi-center classification accuracy (86.1%) was driven by organ and platform confounding. Under controlled, same-cancer same-platform evaluation on Breast Xenium (10x Genomics vs. NCBI GEO repositories), acquisition-domain predictability collapsed to chance:
  $$\text{AUROC} = 0.4900 \quad (95\% \text{ CI: } [0.4422, 0.5364], \text{Balanced Accuracy} = 46.33\%)$$
- **Nomenclature**: Formally designated as **dataset-origin / acquisition-domain prediction** (or acquisition batch), reflecting repository metadata origin rather than independent medical hospital centers.

### 4. Cross-Platform Support Harmonization & Exploratory Variance Partitioning
- **Support Convolution**: Aggregating single-cell Xenium transcripts into $55\,\mu\text{m}$ hexagonal pseudo-spots numerically restored gene PCC from $0.00840$ to $0.01649$, establishing equivalence with the within-Visium reference ($0.01605$; TOST $p < 0.0001$ at margin $\pm 0.003$).
- **Exploratory Linear ANOVA Partitioning (LMG Metric)** with 2,000 bootstrap resamples:
  - **Platform Detection Chemistry (In situ vs. Capture)**: **32.4%** (95% CI: [26.0%, 38.8%])
  - **Physical Support Mismatch (15 $\mu$m vs. 55 $\mu$m)**: **26.8%** (95% CI: [20.9%, 33.3%])
  - **Organ / Biological Divergence**: **25.1%** (95% CI: [19.2%, 31.4%])
  - **Residual / Acquisition Batch Signatures**: **15.7%** (95% CI: [11.0%, 21.2%]) (pooling unmodeled batch signatures and higher-order interactions).

### 5. TOST Equivalence Testing
- At margin $\pm 0.003$: Image-Only vs. Coordinate C1 difference is non-significant ($p_{\text{sup}} = 0.9154$), but equivalence is **NOT** established ($p_{\text{TOST}} = 0.4016$).
- Methodological rule: Non-significant difference does not establish equivalence; deep vision models fail to statistically outperform coordinates while failing equivalence criteria.

### 6. The Proposed VST-Eval Checklist
Rather than claiming a formal societal reporting standard, this benchmark proposes the **VST-Eval Checklist**:
1. *Biological Donor Provenance*: Explicit donor-level tracking and serial section quarantine.
2. *Leakage Ladder Partitioning*: Mandatory patient-held-out (Level 3) and cohort-held-out (Level 4) evaluation.
3. *Spatial Coordinate Nulls*: Obligatory comparison against 2D coordinate baselines (C0–C5).
4. *Variogram-Preserving Spatial Permutations*: Torus or block spatial nulls with explicit test denominators.
5. *Formal Equivalence Testing*: TOST at calibrated biological margins ($\pm 0.003$).
6. *Physical Support Harmonization*: Explicit point-spread convolution matching physical spot PSF before cross-platform inference.
7. *Epistemic Uncertainty Reporting*: Predictive variance calibration and selective prediction curves.

---

## Target Submission Track: Nature Methods (Analysis)
In adherence to *Nature Methods* editorial guidelines, comprehensive empirical evaluations of existing computational methods across multiple datasets and benchmark criteria are categorized as **Analysis** articles (exemplified by recent milestone benchmarks such as the 2025 analysis of 27 single-cell perturbation prediction methods across 29 datasets). This work is structured and locked as an Analysis article providing objective methodological scrutiny.

---

## Citation
```bibtex
@article{vst_eval_integrity_2026,
  title={The Evaluation-Integrity Deficit in Virtual Spatial Transcriptomics: Auditing Data Leakage, Spatial Priors, and Hidden Shortcuts Across Routine Histology},
  author={chessyuan},
  journal={Nature Methods (In Submission, Analysis Track)},
  year={2026}
}
```
