# Soil Microbiome Functional Comparison (16S vs Shotgun)

This repository contains all scripts, data structure, and results used for my MSc thesis:
**Functional Comparison of Soil Microbiomes (Rice vs Shrimp) via 16S and Shotgun Metagenomics**.

Both datasets originate from **[BioProject PRJNA385949 – Soil microbiome functional benchmark](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA385949)**.

---

## 📂 Project Structure

Soil-Microbiome-Pipeline/
├── README.md
├── picrust2_experiment/
│ ├── soil_rep_seqs.fna # ASV sequences (FASTA)
│ ├── soil_otu_table.tsv # ASV abundance table (QIIME2 export)
│ ├── picrust2_out_min10_v2/ # Final corrected PICRUSt2 outputs (KO, EC, MetaCyc, NSTI)
│ ├── metadata.tsv # Example metadata
│ └── picrust2_pipeline.sh # Main PICRUSt2 script (correct FASTA + table inputs)
├── shotgun_rice/
│ ├── SRR5259832_1.fastq.gz # Raw reads (not uploaded; downloadable via SRA)
│ ├── SRR5259832_2.fastq.gz
│ ├── humann_out_SRR5259832/ # HUMAnN output (gene families, KO, RXN, pathways)
│ └── humann_pipeline.sh # HUMAnN run + post-processing
├── scripts/
│ ├── compare_humann_picrust2.py # KO comparison (HUMAnN vs PICRUSt2)
│ └── visualization_pipeline.py # Plotting top EC/RXN functions
└── results/
├── tables/
└── figures/
