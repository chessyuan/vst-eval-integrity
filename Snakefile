# Snakefile for VST-Eval-Integrity
# Unified Reproducible Computational Benchmark Workflow

configfile: "configs/benchmark_config.yaml"

rule all:
    input:
        "donor_audit/DONOR_PROVENANCE_MASTER.tsv",
        "splits/leakage_ladder_metrics.tsv",
        "splits/patient_paired_eor.tsv",
        "coordinate_nulls/coordinate_null_suite.tsv",
        "coordinate_nulls/factorial_image_coord.tsv",
        "spatial_nulls/spatial_autocorrelation_null.tsv",
        "equivalence_tests/tost_equivalence_tests.tsv",
        "equivalence_tests/model_scaling_mixed_effects.tsv",
        "support_harmonization/support_harmonization_ablation.tsv",
        "support_harmonization/variance_decomposition_platform.tsv",
        "uncertainty/table5_uncertainty_failure_detection.tsv",
        "model_adapters/MODEL_STATUS_20260926.tsv",
        "reproduce_figures/fig1_leakage_ladder_eor.png",
        "reproduce_figures/fig2_coordinate_null_tost.png",
        "reproduce_figures/fig3_spatial_permutation_nulls.png",
        "reproduce_figures/fig4_support_harmonization.png",
        "reproduce_figures/fig5_model_scaling_deficit.png",
        "reproduce_figures/fig6_pathway_harmonized_gain.png"

rule donor_audit:
    output:
        "donor_audit/DONOR_PROVENANCE_MASTER.tsv"
    shell:
        "python donor_audit/audit_donor_provenance.py"

rule leakage_ladder:
    input:
        "donor_audit/DONOR_PROVENANCE_MASTER.tsv"
    output:
        "splits/leakage_ladder_metrics.tsv",
        "splits/patient_paired_eor.tsv"
    shell:
        "python splits/leakage_ladder.py"

rule coordinate_nulls:
    input:
        "splits/leakage_ladder_metrics.tsv"
    output:
        "coordinate_nulls/coordinate_null_suite.tsv",
        "coordinate_nulls/factorial_image_coord.tsv"
    shell:
        "python coordinate_nulls/coordinate_models.py"

rule spatial_nulls:
    input:
        "coordinate_nulls/coordinate_null_suite.tsv"
    output:
        "spatial_nulls/spatial_autocorrelation_null.tsv"
    shell:
        "python spatial_nulls/spatial_permutation.py"

rule equivalence_tests:
    input:
        "coordinate_nulls/coordinate_null_suite.tsv"
    output:
        "equivalence_tests/tost_equivalence_tests.tsv",
        "equivalence_tests/model_scaling_mixed_effects.tsv"
    shell:
        "python equivalence_tests/tost_engine.py"

rule support_harmonization:
    output:
        "support_harmonization/support_harmonization_ablation.tsv",
        "support_harmonization/variance_decomposition_platform.tsv"
    shell:
        "python support_harmonization/hexagonal_binning.py"

rule uncertainty_audit:
    output:
        "uncertainty/table5_uncertainty_failure_detection.tsv"
    shell:
        "python uncertainty/selective_prediction.py"

rule model_adapters:
    output:
        "model_adapters/MODEL_STATUS_20260926.tsv"
    shell:
        "python model_adapters/adapters.py"

rule reproduce_figures:
    input:
        "splits/patient_paired_eor.tsv",
        "coordinate_nulls/coordinate_null_suite.tsv",
        "spatial_nulls/spatial_autocorrelation_null.tsv",
        "support_harmonization/support_harmonization_ablation.tsv",
        "equivalence_tests/model_scaling_mixed_effects.tsv"
    output:
        "reproduce_figures/fig1_leakage_ladder_eor.png",
        "reproduce_figures/fig2_coordinate_null_tost.png",
        "reproduce_figures/fig3_spatial_permutation_nulls.png",
        "reproduce_figures/fig4_support_harmonization.png",
        "reproduce_figures/fig5_model_scaling_deficit.png",
        "reproduce_figures/fig6_pathway_harmonized_gain.png"
    shell:
        "python reproduce_figures/generate_benchmark_figures.py"

rule benchmark:
    input:
        rules.all.input
