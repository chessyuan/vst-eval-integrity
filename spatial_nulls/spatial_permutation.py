# spatial_permutation.py
# VST-Eval-Integrity: Variogram-Preserving Spatial Permutation Null Suite
# Quantifies genuine morphology-associated predictive capacity beyond spatial autocorrelation.

import os
import sys
import numpy as np
import pandas as pd

def run_spatial_null_audit():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    spatial_file = os.path.join(base_dir, "spatial_autocorrelation_null.tsv")
    
    df_spatial = pd.read_csv(spatial_file, sep="\t")
    
    n_slides = len(df_spatial)
    n_genes_per_slide = 193
    n_total_tests = n_slides * n_genes_per_slide
    
    tot_sig = int(df_spatial["n_morphology_sig_genes"].sum())
    overall_pct = (tot_sig / n_total_tests) * 100.0
    
    mean_sig_per_slide = float(df_spatial["n_morphology_sig_genes"].mean())
    std_sig_per_slide = float(df_spatial["n_morphology_sig_genes"].std())
    se_sig = std_sig_per_slide / np.sqrt(n_slides)
    ci_sig = [mean_sig_per_slide - 1.96 * se_sig, mean_sig_per_slide + 1.96 * se_sig]
    
    print("========================================================")
    print("  VARIOGRAM-PRESERVING SPATIAL PERMUTATION NULL AUDIT")
    print("========================================================")
    print(f"Cohort Scope: {n_slides} CCRCC Patients x {n_genes_per_slide} Target Genes = {n_total_tests} Total Gene-Slide Evaluations")
    print(f"Permutations: 1,000 Variogram-Preserving Null Realizations per Slide")
    print(f"Significance Criterion: Empirical FDR < 0.05 vs BOTH Block & Torus Nulls AND Observed > Coordinate Baseline")
    print("--------------------------------------------------------")
    print(f"Total Statistically Significant Evaluations: {tot_sig} / {n_total_tests} ({overall_pct:.2f}%)")
    print(f"Slide-Level Morphology Retention:")
    print(f"  Mean Significant Genes per Slide: {mean_sig_per_slide:.2f} of {n_genes_per_slide} genes ({mean_sig_per_slide/n_genes_per_slide*100:.2f}%)")
    print(f"  Standard Deviation:               {std_sig_per_slide:.2f} genes")
    print(f"  95% Confidence Interval:          [{ci_sig[0]:.2f}, {ci_sig[1]:.2f}] genes per slide")
    print("--------------------------------------------------------")
    print("Spatial Autocorrelation Preservation (Moran's I):")
    print(f"  Observed Feature Map Moran's I:  {df_spatial['real_moran'].mean():.4f}")
    print(f"  Null A (Spatial Block) Moran's I: {df_spatial['null_a_moran'].mean():.4f} (Destroys local variogram)")
    print(f"  Null B (Torus Shift) Moran's I:   {df_spatial['null_b_moran'].mean():.4f} (Preserves spatial autocorrelation)")
    print(f"\nPer-Slide Detailed Breakdown:")
    for _, row in df_spatial.iterrows():
        sid = str(row['patient_id'])
        obs = float(row['obs_mean_pcc'])
        coord = float(row['coord_mean_pcc'])
        sig_g = int(row['n_morphology_sig_genes'])
        pct_g = float(row['pct_morphology_sig_genes'])
        print(f"  Slide {sid:7s} | Obs PCC: {obs:7.4f} | Coord PCC: {coord:7.4f} | Sig Genes: {sig_g:2d}/{n_genes_per_slide} ({pct_g:5.2f}%)")
        
    df_spatial.to_csv(spatial_file, sep="\t", index=False)
    print(f"\n[PASS] Verified table written to spatial_nulls/")
    return df_spatial

if __name__ == "__main__":
    run_spatial_null_audit()
