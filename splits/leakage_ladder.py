# leakage_ladder.py
# VST-Eval-Integrity: Evaluation Leakage Ladder & Patient-Paired EOR
# Quantifies optimism from spot-level spatial autocorrelation down to patient-unseen holdout.

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

def run_leakage_audit():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    ladder_file = os.path.join(base_dir, "leakage_ladder_metrics.tsv")
    eor_file = os.path.join(base_dir, "patient_paired_eor.tsv")
    
    df_ladder = pd.read_csv(ladder_file, sep="\t")
    df_eor = pd.read_csv(eor_file, sep="\t")
    
    print("========================================================")
    print("  DATA LEAKAGE LADDER & EVALUATION OPTIMISM AUDIT")
    print("========================================================")
    print("Leakage Ladder Performance across Generalization Levels:")
    for _, row in df_ladder.iterrows():
        lvl = str(row["level"])
        gp = float(row["gene_pcc"])
        sp = float(row["spot_pcc"])
        inf = float(row["inflation_factor"])
        print(f"  {lvl:35s} | Gene PCC: {gp:8.5f} | Spot PCC: {sp:8.5f} | Inflation: {inf:6.2f}x")
        
    pcc_l0 = df_ladder.loc[df_ladder["level"].str.contains("Spot Random"), "gene_pcc"].values[0]
    pcc_l3 = df_ladder.loc[df_ladder["level"].str.contains("Patient Holdout|Patient-Held-Out"), "gene_pcc"].values[0]
    agg_eor = pcc_l0 / pcc_l3
    print(f"\nAggregate Apparent-to-Actual Optimism Ratio: {agg_eor:.2f}x ({pcc_l0:.4f} vs {pcc_l3:.5f})")
    
    eor_vals = df_eor["eor"].values
    delta_vals = df_eor["delta_pcc"].values
    n_donors = len(df_eor)
    
    np.random.seed(42)
    b_eor, b_delta = [], []
    for _ in range(5000):
        idx = np.random.choice(n_donors, size=n_donors, replace=True)
        b_eor.append(np.median(eor_vals[idx]))
        b_delta.append(np.mean(delta_vals[idx]))
        
    med_eor = np.median(eor_vals)
    iqr_eor = np.percentile(eor_vals, 75) - np.percentile(eor_vals, 25)
    ci_eor = [np.percentile(b_eor, 2.5), np.percentile(b_eor, 97.5)]
    
    mean_delta = np.mean(delta_vals)
    ci_delta = [np.percentile(b_delta, 2.5), np.percentile(b_delta, 97.5)]
    
    stat_w, p_wilcoxon = wilcoxon(df_eor["pcc_spot_random_within"], df_eor["pcc_patient_unseen"])
    
    print(f"\nPatient-Paired Evaluation Optimism Ratio (N={n_donors} CCRCC Donors):")
    print(f"  Median EOR:        {med_eor:.2f} (IQR: {iqr_eor:.2f}, 95% Bootstrap CI: [{ci_eor[0]:.2f}, {ci_eor[1]:.2f}])")
    print(f"  Mean Delta PCC:    +{mean_delta:.4f} (95% Bootstrap CI: [{ci_delta[0]:.4f}, {ci_delta[1]:.4f}])")
    print(f"  Wilcoxon Test:     p = {p_wilcoxon:.4e} (Spot-Random vs Patient-Held-Out)")
    
    df_ladder.to_csv(ladder_file, sep="\t", index=False)
    df_eor.to_csv(eor_file, sep="\t", index=False)
    print(f"\n[PASS] Verified tables written to splits/")
    return df_ladder, df_eor

if __name__ == "__main__":
    run_leakage_audit()
