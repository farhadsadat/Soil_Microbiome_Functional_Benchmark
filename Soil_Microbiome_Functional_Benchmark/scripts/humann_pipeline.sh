#!/bin/bash
# HUMAnN full pipeline

# Preprocess (fastp, reformat, etc.)
fastp -i SRR5259832_1.fastq.gz -I SRR5259832_2.fastq.gz \
      -o SRR5259832_1.clean.fq.gz -O SRR5259832_2.clean.fq.gz \
      --detect_adapter_for_pe --html fastp_report.html

reformat.sh in1=SRR5259832_1.clean.fq.gz in2=SRR5259832_2.clean.fq.gz \
            out=SRR5259832_merged.clean.fastq minlen=30

# Taxonomic profiling
metaphlan SRR5259832_merged.clean.fastq \
    --input_type fastq --nproc 8 \
    -o SRR5259832_metaphlan_bugs_list.tsv

# HUMAnN execution
humann --input SRR5259832_merged.clean.fastq \
      --output ../humann_out_SRR5259832 \
      --threads 8 \
      --nucleotide-database ~/humann_db/chocophlan \
      --protein-database ~/humann_db/uniref90 \
      --taxonomic-profile SRR5259832_metaphlan_bugs_list.tsv

# Normalization
humann_renorm_table --input ../humann_out_SRR5259832/SRR5259832_merged.clean_genefamilies.tsv \
                    --output ../humann_out_SRR5259832/SRR5259832_genefamilies_relab.tsv \
                    --units cpm

# Regroup to KO (custom map)
humann_regroup_table --input ../humann_out_SRR5259832/SRR5259832_genefamilies_relab.tsv \
                     --custom ~/humann_db/utility_mapping/map_uniref90_ko_clean.tsv.gz \
                     --output ../humann_out_SRR5259832/SRR5259832_KO_final.tsv \
                     --function sum

humann_split_stratified_table --input ../humann_out_SRR5259832/SRR5259832_KO_final.tsv \
                              --output ../humann_out_SRR5259832/SRR5259832_KO_final

humann_renorm_table --input ../humann_out_SRR5259832/SRR5259832_KO_final/SRR5259832_KO_final_unstratified.tsv \
                    --output ../humann_out_SRR5259832/SRR5259832_KO_unstrat_cpm.tsv \
                    --units cpm
