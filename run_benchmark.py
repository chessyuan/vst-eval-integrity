# -*- coding: utf-8 -*-
import os, sys, time, subprocess

PIPELINE_STAGES = [
    ('1. Biological Donor Provenance Audit', 'donor_audit/audit_donor_provenance.py'),
    ('2. 7-Level Leakage Ladder & Patient-Paired EOR', 'splits/leakage_ladder.py'),
    ('3. Coordinate-Only Null Suite (C0-C5) & Factorial Analysis', 'coordinate_nulls/coordinate_models.py'),
    ('4. Variogram-Preserving Spatial Permutation Nulls', 'spatial_nulls/spatial_permutation.py'),
    ('5. TOST Equivalence Testing & Model Scaling Mixed Effects', 'equivalence_tests/tost_engine.py'),
    ('6. Support Harmonization & Exploratory ANOVA Partitioning', 'support_harmonization/hexagonal_binning.py'),
    ('7. Epistemic Uncertainty & Selective Prediction Retention', 'uncertainty/selective_prediction.py'),
    ('8. Unified Model Inference Adapters Audit', 'model_adapters/adapters.py'),
    ('9. Reproduce All 6 Publication Benchmark Figures', 'reproduce_figures/generate_benchmark_figures.py')
]

def main():
    print('=' * 75)
    print('  VST-EVAL-INTEGRITY: REPRODUCIBLE COMPUTATIONAL BENCHMARK HARNESS')
    print('  Target Venue: Nature Methods (Analysis Track)')
    print('  Maintainer: Repository maintainer')
    print('=' * 75)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    total_start = time.time()
    for idx, (title, script_rel) in enumerate(PIPELINE_STAGES, 1):
        script_path = os.path.join(base_dir, script_rel)
        print(f'\n[{idx}/{len(PIPELINE_STAGES)}] Running: {title}')
        print(f'    Script: {script_rel}')
        stage_start = time.time()
        res = subprocess.run([sys.executable, script_path], cwd=base_dir, capture_output=True, text=True)
        stage_elapsed = time.time() - stage_start
        if res.returncode == 0:
            print(f'    Status: [PASS] ({stage_elapsed:.2f} seconds)')
        else:
            print(f'    Status: [FAIL] (code {res.returncode})')
            print('    STDERR:', res.stderr[:500])
            sys.exit(1)
    total_elapsed = time.time() - total_start
    print('\n' + '=' * 75)
    print(f'  [PASS] ALL {len(PIPELINE_STAGES)} BENCHMARK STAGES COMPLETED SUCCESSFULLY ({total_elapsed:.2f} seconds)')
    print('=' * 75)

if __name__ == '__main__':
    main()
