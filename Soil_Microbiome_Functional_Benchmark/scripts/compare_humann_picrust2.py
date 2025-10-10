#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

ap = argparse.ArgumentParser()
ap.add_argument("--humann", required=True, help="HUMAnN KO unstrat (e.g., *_KO_final_unstratified.tsv or *_KO_cpm.tsv)")
ap.add_argument("--picrust", default="picrust2_experiment/KO_unstrat_annot.tsv.gz",
                help="PICRUSt2 KO table (columns: function, description, soil_rice, soil_shrimp)")
ap.add_argument("--picrust-sample", default="soil_rice", help="Which PICRUSt2 sample column to compare")
ap.add_argument("--out", default="results/tables/KO_compare.tsv", help="Output merged table")
ap.add_argument("--png", default="results/figures/KO_scatter.png", help="Output scatter plot")
args = ap.parse_args()

def read_table(path, **kw):
    comp = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, sep="\t", compression=comp, **kw)

# HUMAnN: two columns if unstratified; drop comments
h = read_table(args.humann, header=None, comment="#")
# Try to detect header vs no-header
if h.shape[1] >= 2 and not h.iloc[0,0].startswith("K"):
    # has header; rename if needed
    h = read_table(args.humann, sep="\t", header=0)
    # first col name may be like "# Gene Family"
    first = h.columns[0]
    h = h[[first, h.columns[1]]]
    h.columns = ["KO", "HUMANN"]
else:
    h = h.iloc[:, :2]
    h.columns = ["KO", "HUMANN"]

h = h[h["KO"].astype(str).str.match("^K")]

# PICRUSt2: has columns: function, description, soil_rice, soil_shrimp
p = read_table(args.picrust)
p.columns = [c.strip() for c in p.columns]
p = p.rename(columns={"function":"KO"})
p["KO"] = p["KO"].str.replace("ko:","", regex=False)
if args.picrust_sample not in p.columns:
    raise SystemExit(f"Column {args.picrust_sample} not in {list(p.columns)}")
p = p[["KO", args.picrust_sample]].rename(columns={args.picrust_sample: "PICRUSt2"})

# Merge
merged = pd.merge(h, p, on="KO", how="inner")
merged = merged.replace([np.inf, -np.inf], np.nan).dropna()

print(f"Merged {len(merged)} shared KOs (HUMAnN:{len(h)}, PICRUSt2:{len(p)})")
if len(merged) >= 2:
    s = spearmanr(merged["HUMANN"], merged["PICRUSt2"])
    r = pearsonr(merged["HUMANN"], merged["PICRUSt2"])
    print(f"Spearman r={s.statistic:.3f} (p={s.pvalue:.2e}); Pearson r={r.statistic:.3f} (p={r.pvalue:.2e})")

# Plot
eps = 1e-9
plt.figure(figsize=(6,5))
plt.scatter(np.log10(merged["PICRUSt2"]+eps), np.log10(merged["HUMANN"]+eps), alpha=0.7, s=22)
plt.xlabel(f"log10 PICRUSt2 ({args.picrust_sample})")
plt.ylabel("log10 HUMAnN")
plt.title("KO correlation: HUMAnN vs PICRUSt2")
plt.tight_layout()
os.makedirs("results/figures", exist_ok=True)
plt.savefig(args.png, dpi=300)
print(f"Saved {args.png}")

os.makedirs("results/tables", exist_ok=True)
merged.to_csv(args.out, sep="\t", index=False)
print(f"Wrote {args.out}")
