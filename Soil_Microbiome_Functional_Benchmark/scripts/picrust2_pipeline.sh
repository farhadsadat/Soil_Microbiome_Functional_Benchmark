#!/usr/bin/env bash
set -euo pipefail

# === PICRUSt2 PIPELINE (16S soil microbiome) ===
# Folder structure (inside project):
# picrust2_experiment/
# ├── soil_rep_seqs.fna
# ├── soil_otu_table.tsv
# └── picrust2_out_min10_v2/

# === Define inputs and outputs ===
INPUT_DIR="picrust2_experiment"
ASV_FASTA="${INPUT_DIR}/soil_rep_seqs.fna"          # FASTA of ASV sequences
ASV_TABLE="${INPUT_DIR}/soil_otu_table.tsv"         # ASV abundance table
OUT_DIR="${INPUT_DIR}/picrust2_out_min10_v2"        # Output directory
THREADS="${THREADS:-8}"

mkdir -p "${OUT_DIR}"

echo "------------------------------------------------------"
echo "[PICRUSt2] Starting functional prediction"
echo " Input FASTA : ${ASV_FASTA}"
echo " Input Table : ${ASV_TABLE}"
echo " Output Dir  : ${OUT_DIR}"
echo " Threads     : ${THREADS}"
echo "------------------------------------------------------"

# === 1. Run PICRUSt2 core pipeline ===
picrust2_pipeline.py \
  -s "${ASV_FASTA}" \
  -i "${ASV_TABLE}" \
  -o "${OUT_DIR}" \
  -p "${THREADS}"

echo
echo "[PICRUSt2] Run complete. Generating NSTI summary..."

# === 2. Quick NSTI summary ===
NSTI_FILE="${OUT_DIR}/marker_nsti/nsti.tsv"
if [[ -f "${NSTI_FILE}" ]]; then
  awk 'NR>1{print $2}' "${NSTI_FILE}" | \
  awk '{sum+=$1;a[NR]=$1}END{n=NR;asort(a);printf("NSTI  n=%d  mean=%.4f  median=%.4f  p90=%.4f  p95=%.4f\n",n,sum/n,a[int(n/2)],a[int(0.9*n)],a[int(0.95*n)])}'
else
  echo "[!] NSTI file not found (${NSTI_FILE})"
fi

echo
echo "[✓] PICRUSt2 complete — results saved in:"
echo "    ${OUT_DIR}"
echo "------------------------------------------------------"
echo "Expected output files:"
echo " • KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
echo " • EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
echo " • pathways_out/pathway_abundance.tsv.gz"
echo " • marker_nsti/nsti.tsv"
echo "------------------------------------------------------"
