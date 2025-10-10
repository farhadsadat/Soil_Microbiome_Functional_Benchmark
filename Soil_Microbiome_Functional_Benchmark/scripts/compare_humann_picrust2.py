import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import numpy as np

humann = "../humann_out_SRR5259832/HUMANN_KO_manual_cpm.tsv"
picrust = "../KO_unstrat_annot.tsv.gz"

h = pd.read_csv(humann, sep="\\t", comment="#", header=None, names=["KO","HUMAnN"])
h = h[h["KO"].str.match("^K")]

p = pd.read_csv(picrust, sep="\\t")
p = p.rename(columns={"function":"KO"})
p["KO"] = p["KO"].str.replace("ko:","", regex=False)
p = p[["KO","soil_rice"]].rename(columns={"soil_rice":"PICRUSt2"})

merged = pd.merge(h, p, on="KO", how="inner")

if len(merged) > 1:
    s = spearmanr(merged["HUMAnN"], merged["PICRUSt2"])
    pcc = pearsonr(merged["HUMAnN"], merged["PICRUSt2"])
    print("Spearman:", s.statistic, "p=", s.pvalue)
    print("Pearson:", pcc.statistic, "p=", pcc.pvalue)

plt.figure(figsize=(6,5))
plt.scatter(np.log10(merged["PICRUSt2"] + 1e-9), np.log10(merged["HUMAnN"] + 1e-9), alpha=0.7)
plt.xlabel("log10 PICRUSt2")
plt.ylabel("log10 HUMAnN")
plt.title("KO correlation: HUMAnN vs PICRUSt2")
plt.tight_layout()
plt.savefig("KO_scatter_HUMAnN_vs_PICRUSt2.png", dpi=300)
