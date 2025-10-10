# Soil_Microbiome_Functional_Benchmark



**Author:** Farhad Sadat
**Affiliation:** University of Padova 
**Project Type:** Master’s Thesis – Comparative Functional Profiling  
**Tools:** HUMAnN 3.9, PICRUSt2, MetaPhlAn, Python (pandas/matplotlib)  
**Objective:** Benchmark shotgun metagenomic functional profiling (HUMAnN) against 16S-based predictive profiling (PICRUSt2) using soil samples from rice and shrimp farming environments.

---

## 📖 Overview

This repository contains all scripts, tables, and figures used to compare microbial functional inference methods in soil metagenomes.  
The goal was to evaluate how well *inferred* 16S-based functions (PICRUSt2) correlate with *directly observed* shotgun-based functions (HUMAnN3).

---

## 📂 Repository Structure


---

## 🧪 Datasets

### 🧫 Shotgun Metagenome
- **SRA Run:** [SRR5259832](https://www.ncbi.nlm.nih.gov/sra/SRR5259832)  
- **BioProject:** [PRJNA385949](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA385949)  
- **Description:** Soil metagenomic reads from rice paddy fields in Asia.  
- **Used for:** Direct functional annotation via **HUMAnN 3.9**.

### 🧬 16S rRNA Dataset
- **SRA Run:** [PRJNA639700](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639700)  
- **Description:** 16S rRNA gene amplicons from rice and shrimp soil microbiomes.  
- **Used for:** Functional prediction via **PICRUSt2**.

---

## 🧠 Methodology Summary

| Step | Tool / Method | Description |
|------|----------------|-------------|
| 1 | **Quality control** | Performed using `fastp` for adapter trimming and quality filtering. |
| 2 | **Taxonomic profiling** | Run `MetaPhlAn` on merged reads to estimate community composition. |
| 3 | **Functional profiling** | HUMAnN pipeline to map reads to UniRef90 → KO, EC, and MetaCyc reactions. |
| 4 | **16S prediction** | PICRUSt2 used on ASV/OTU tables to infer KO and pathways. |
| 5 | **Benchmarking** | Compared HUMAnN vs PICRUSt2 KO abundances (Spearman, Pearson). |
| 6 | **Visualization** | Generated bar charts, heatmaps, and EC class composition using Python. |

---

## 📜 Scripts

### 🧬 `humann_pipeline.sh`
Runs the complete HUMAnN workflow from raw shotgun reads to KO tables.  
Includes:
- Quality control (`fastp`, `reformat.sh`)
- Taxonomic profiling (`MetaPhlAn`)
- Functional profiling (`HUMAnN`)
- Regrouping (UniRef → KO)
- Normalization (CPM)
- Output splitting and renormalization

### 🧫 `picrust2_pipeline.sh`
Executes the PICRUSt2 pipeline starting from 16S OTU/ASV tables.  
Outputs unstratified KO and pathway tables.

### ⚖️ `compare_humann_picrust2.py`
Merges HUMAnN and PICRUSt2 KO tables to calculate:
- Spearman and Pearson correlation between methods
- Scatter plots (`KO_scatter_HUMAnN_vs_PICRUSt2.png`)
- KO overlap statistics (`KO_overlap_HUMAnN_vs_PICRUSt2.tsv`)

### 📊 `visualization_pipeline.py`
Generates all HUMAnN-based figures:
- `EC_top20.png`, `RXN_top20.png`
- `EC_heatmap.png`, `RXN_heatmap.png`
- `EC_class_stacked.png`
  
Usage:
```bash
python scripts/visualization_pipeline.py \
  --ec  ~/shotgun_rice/humann_out_SRR5259832/EC_manual_cpm_named.tsv \
  --rxn ~/shotgun_rice/humann_out_SRR5259832/SRR5259832_RXN_cpm_named.tsv \
  --topn 20 \
  --outdir figures
