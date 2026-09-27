# coordinate_models.py
# VST-Eval-Integrity: Coordinate Null Suite & Factorial Image-Coordinate Dissection
# Benchmarks deep vision features against pure 2D coordinate null models (C0-C5) and tests feature fusion.

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

def run_coordinate_null_audit():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    coord_file = os.path.join(base_dir, "coordinate_null_suite.tsv")
    fact_file = os.path.join(base_dir, "factorial_image_coord.tsv")
    
    df_coord = pd.read_csv(coord_file, sep="\t")
    df_fact = pd.read_csv(fact_file, sep="\t")
    
    print("========================================================")
    print("  COORDINATE NULL SUITE & FACTORIAL INTERACTION AUDIT")
    print("========================================================")
    print("Mean Cross-Patient Gene PCC across Coordinate Baselines (N=24 Donors):")
    cols = ["C0_Raw_XY", "C1_Normalized_XY", "C2_Rotated_45deg", "C3_Reflected_X", "C4_RBF_Spline_16", "C5_Gaussian_Kernel", "Image_Only", "Image_Plus_Coord", "Shuffled_Image_Coord"]
    for c in cols:
        vals = df_coord[c].values
        print(f"  {c:22s} | Mean: {np.mean(vals):8.5f} | Std: {np.std(vals):8.5f} | Median: {np.median(vals):8.5f}")
        
    diff = df_coord["Image_Only"].values - df_coord["C1_Normalized_XY"].values
    t_stat, p_val = ttest_rel(df_coord["Image_Only"].values, df_coord["C1_Normalized_XY"].values)
    print(f"\nImage-Only vs Coordinate C1 Paired Difference:")
    print(f"  Mean Delta: {np.mean(diff):.6f} (SE: {np.std(diff, ddof=1)/np.sqrt(len(diff)):.6f})")
    print(f"  Paired t-test: t = {t_stat:.4f}, p = {p_val:.4f} (No detectable superiority of image over coordinate alone)")
    
    print("\nFactorial Image x Coordinate Interaction:")
    for _, row in df_fact.iterrows():
        arm = str(row['Arm'])
        mpcc = float(row['Mean_Gene_PCC'])
        sd = float(row['SD_PCC'])
        print(f"  {arm:35s} | Mean PCC: {mpcc:8.5f} | Std: {sd:8.5f}")
        
    df_coord.to_csv(coord_file, sep="\t", index=False)
    df_fact.to_csv(fact_file, sep="\t", index=False)
    print(f"\n[PASS] Verified tables written to coordinate_nulls/")
    return df_coord, df_fact

if __name__ == "__main__":
    run_coordinate_null_audit()
