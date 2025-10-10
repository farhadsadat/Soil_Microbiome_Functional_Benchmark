#!/bin/bash
# PICRUSt2 pipeline (16S → KO/Pathway predictions)

# Input: soil_otu_table.tsv or 16S FastQ (preprocessed)
# Example using PICRUSt2:
picrust2_pipeline.py -s soil_otu_table.tsv -o picrust2_out \
    --output-kegg KO_unstrat_annot.tsv \
    --output-pathway PATH_unstrat_annot.tsv \
    --threads 8 \
    --verbose

# (Optional: gzip the ko table)
gzip -f KO_unstrat_annot.tsv
