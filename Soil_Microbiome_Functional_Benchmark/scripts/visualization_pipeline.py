#!/usr/bin/env python3
"""
Visualization pipeline for HUMAnN outputs (EC & MetaCyc RXNs)

Generates:
- EC_top20.png               (bar chart of top 20 ECs by CPM)
- RXN_top20.png              (bar chart of top 20 RXNs by CPM)
- EC_heatmap.png             (heatmap of top 20 ECs, row z-score of log10(CPM))
- RXN_heatmap.png            (heatmap of top 20 RXNs, row z-score of log10(CPM))
- EC_class_stacked.png       (stacked bar of EC classes 1..6 + other)

Inputs (defaults can be overridden with CLI flags):
- EC table (CPM) with optional header and name column:
    ~/shotgun_rice/humann_out_SRR5259832/EC_manual_cpm_named.tsv
    Expected columns (minimum): EC_ID, CPM [, Name]
- RXN table (CPM) with optional header and name column:
    ~/shotgun_rice/humann_out_SRR5259832/SRR5259832_RXN_cpm_named.tsv
    Expected columns (minimum): RXN_ID, CPM [, Name]

Author: Farhad Sadat
"""

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------- helpers -------------------------- #

def _read_two_cols(path: str, id_idx=0, val_idx=1):
    """
    Read a TSV that may contain a commented header line.
    Returns a DataFrame with columns ['ID','Abundance'].
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input not found: {path}")

    # Try with header=None then fallback to header=0
    try:
        df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
        # If file actually has a header row, detect non-numeric abundance and reload with header=0
        if df.shape[1] > max(id_idx, val_idx):
            try:
                _ = pd.to_numeric(df.iloc[1, val_idx])
            except Exception:
                df = pd.read_csv(path, sep="\t", comment="#", header=0, dtype=str)
        else:
            df = pd.read_csv(path, sep="\t", comment="#", header=0, dtype=str)
    except Exception:
        df = pd.read_csv(path, sep="\t", comment="#", header=0, dtype=str)

    # Select columns safely
    if df.shape[1] <= max(id_idx, val_idx):
        raise ValueError(f"File {path} does not have the expected number of columns.")

    sub = df.iloc[:, [id_idx, val_idx]].copy()
    sub.columns = ["ID", "Abundance"]
    # coerce abundance to numeric
    sub["Abundance"] = pd.to_numeric(sub["Abundance"], errors="coerce")
    sub = sub.dropna(subset=["ID", "Abundance"])
    # Keep positive values
    sub = sub[sub["Abundance"] > 0]
    return sub


def _topbar(df: pd.DataFrame, title: str, outpng: str, n: int = 20, width: float = 8.0):
    """
    Barh of top n rows by 'Abundance'
    """
    d = df.sort_values("Abundance", ascending=False).head(n)
    plt.figure(figsize=(width, max(4, 0.35 * len(d))))
    plt.barh(d["ID"][::-1], d["Abundance"][::-1])
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()
    print(f"Saved {outpng}")


def _zscore_rows(vec: np.ndarray) -> np.ndarray:
    mu = np.nanmean(vec)
    sd = np.nanstd(vec)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(vec)
    return (vec - mu) / sd


def _heatmap_top(df: pd.DataFrame, title: str, outpng: str, n: int = 20):
    d = df.sort_values("Abundance", ascending=False).head(n).copy()
    # log10 transform to compress dynamic range
    eps = 1e-12
    logv = np.log10(d["Abundance"].values + eps)
    z = _zscore_rows(logv).reshape(-1, 1)

    plt.figure(figsize=(6, max(4, 0.35 * len(d))))
    plt.imshow(z, aspect="auto")
    plt.colorbar(label="Row z-score (log10 CPM)")
    plt.yticks(range(len(d)), d["ID"])
    plt.xticks([0], ["Sample"])
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()
    print(f"Saved {outpng}")


def _ec_class(ec_id: str) -> str:
    m = re.match(r"^([1-6])\.", str(ec_id))
    return m.group(1) if m else "other"


def _stacked_ec_classes(ec_df: pd.DataFrame, outpng: str, title: str = "EC top-level class composition"):
    tmp = ec_df.copy()
    tmp["class"] = tmp["ID"].map(_ec_class)
    cls = tmp.groupby("class", as_index=False)["Abundance"].sum()
    # Order classes 1..6 then other if present
    order = [c for c in ["1", "2", "3", "4", "5", "6", "other"] if c in set(cls["class"])]
    cls = cls.set_index("class").loc[order].reset_index()

    plt.figure(figsize=(6, 4))
    bottom = 0.0
    for _, row in cls.iterrows():
        plt.bar([0], [row["Abundance"]], bottom=bottom, label=row["class"])
        bottom += row["Abundance"]
    plt.xticks([0], ["Sample"])
    plt.ylabel("CPM")
    plt.title(title)
    plt.legend(title="EC class", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()
    print(f"Saved {outpng}")

# -------------------------- main -------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Generate HUMAnN visualization figures from EC and RXN CPM tables.")
    ap.add_argument("--ec",  default=os.path.expanduser("~/shotgun_rice/humann_out_SRR5259832/EC_manual_cpm_named.tsv"),
                    help="Path to EC CPM table (default: ~/shotgun_rice/humann_out_SRR5259832/EC_manual_cpm_named.tsv)")
    ap.add_argument("--rxn", default=os.path.expanduser("~/shotgun_rice/humann_out_SRR5259832/SRR5259832_RXN_cpm_named.tsv"),
                    help="Path to RXN CPM table (default: ~/shotgun_rice/humann_out_SRR5259832/SRR5259832_RXN_cpm_named.tsv)")
    ap.add_argument("--topn", type=int, default=20, help="Top N features for bar/heatmaps (default: 20)")
    ap.add_argument("--outdir", default=".", help="Output directory for PNGs (default: current dir)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read
    ec_df  = _read_two_cols(args.ec,  id_idx=0, val_idx=1)
    rxn_df = _read_two_cols(args.rxn, id_idx=0, val_idx=1)

    # Bar charts
    _topbar(ec_df,  f"Top {args.topn} EC numbers (CPM)",   os.path.join(args.outdir, "EC_top20.png"),  n=args.topn)
    _topbar(rxn_df, f"Top {args.topn} MetaCyc reactions (CPM)", os.path.join(args.outdir, "RXN_top20.png"), n=args.topn)

    # Heatmaps
    _heatmap_top(ec_df,  f"Top {args.topn} ECs (z-score of log10 CPM)",   os.path.join(args.outdir, "EC_heatmap.png"),  n=args.topn)
    _heatmap_top(rxn_df, f"Top {args.topn} MetaCyc RXNs (z-score of log10 CPM)", os.path.join(args.outdir, "RXN_heatmap.png"), n=args.topn)

    # Stacked EC classes
    _stacked_ec_classes(ec_df, os.path.join(args.outdir, "EC_class_stacked.png"))

if __name__ == "__main__":
    main()
