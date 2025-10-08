import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# ---------- PATHS ----------
humann = "/Users/farhadsadat/shotgun_rice/humann_out_SRR5259832/HUMANN_KO_manual_cpm.tsv"
picrust = "/Users/farhadsadat/picrust2_experiment/KO_unstrat_annot.tsv.gz"

# ---------- LOAD ----------
# HUMAnN KO table
h = pd.read_csv(humann, sep="\t", comment="#", names=["KO","HUMAnN"])
h = h[h["KO"].str.match("^K")]

# PICRUSt2 KO table (has multiple samples, keep soil_rice only)
p = pd.read_csv(picrust, sep="\t")
p.columns = [c.strip() for c in p.columns]
p = p.rename(columns={"function":"KO"})
p["KO"] = p["KO"].str.replace("ko:","",regex=False)
p = p[["KO","soil_rice"]]        # use one sample (you can change to soil_shrimp later)
p = p.rename(columns={"soil_rice":"PICRUSt2"})

# ---------- MERGE ----------
merged = pd.merge(h, p, on="KO", how="inner")
print(f"Merged {len(merged)} shared KOs out of HUMAnN:{len(h)} / PICRUSt2:{len(p)}")

# ---------- CORRELATIONS ----------
if len(merged) > 1:
    s = spearmanr(merged["HUMAnN"], merged["PICRUSt2"])
    pcc = pearsonr(merged["HUMAnN"], merged["PICRUSt2"])
    print(f"Spearman r={s.statistic:.3f} (p={s.pvalue:.2e})")
    print(f"Pearson  r={pcc.statistic:.3f} (p={pcc.pvalue:.2e})")

# ---------- SCATTER PLOT ----------
plt.figure(figsize=(6,5))
plt.scatter(np.log10(merged["PICRUSt2"]+1e-9),
            np.log10(merged["HUMAnN"]+1e-9),
            alpha=0.7, s=25)
plt.xlabel("log10(PICRUSt2 abundance + 1e-9)")
plt.ylabel("log10(HUMAnN abundance + 1e-9)")
plt.title("KO abundance correlation: HUMAnN vs PICRUSt2 (soil_rice)")
plt.tight_layout()
plt.savefig("KO_scatter_HUMAnN_vs_PICRUSt2.png", dpi=300)
plt.close()

# ---------- TOP-50 OVERLAP ----------
topH = set(h.sort_values("HUMAnN", ascending=False).head(50)["KO"])
topP = set(p.sort_values("PICRUSt2", ascending=False).head(50)["KO"])
overlap = len(topH & topP)
print(f"Top-50 KO overlap: {overlap}/50")

# ---------- SAVE ----------
merged.to_csv("KO_overlap_HUMAnN_vs_PICRUSt2.tsv", sep="\t", index=False)
