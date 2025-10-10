#!/usr/bin/env bash
set -euo pipefail

# conda activate humann39

IN_DIR="shotgun_rice"
SRA="SRR5259832"
OUT_DIR="${IN_DIR}/humann_out_${SRA}"
THREADS="${THREADS:-8}"

# Set your DB paths in the shell before running, e.g.:
# export HUMANN_CHOCOPHLAN=~/humann_db/chocophlan/chocophlan
# export HUMANN_UNIREF=~/humann_db/uniref90/uniref

R1="${IN_DIR}/${SRA}_1.fastq.gz"
R2="${IN_DIR}/${SRA}_2.fastq.gz"
C1="${IN_DIR}/${SRA}_1.clean.fq.gz"
C2="${IN_DIR}/${SRA}_2.clean.fq.gz"
MERGED="${IN_DIR}/${SRA}_merged.clean.fastq"
META="${OUT_DIR}/${SRA}_metaphlan_bugs_list.tsv"

mkdir -p "${IN_DIR}" "${OUT_DIR}"

# 1) QC and merge (keeps your prior layout)
if [[ ! -f "${MERGED}" ]]; then
  echo "[HUMAnN] fastp…"
  fastp -i "${R1}" -I "${R2}" -o "${C1}" -O "${C2}" \
        -h "${IN_DIR}/${SRA}_fastp.html" -j "${IN_DIR}/${SRA}_fastp.json"
  # merge by concatenation (HUMAnN can read single FASTQ)
  zcat "${C1}" "${C2}" > "${MERGED}"
fi

# 2) MetaPhlAn profile (recommended for HUMAnN prescreen)
if [[ ! -f "${META}" ]]; then
  echo "[HUMAnN] MetaPhlAn…"
  metaphlan "${MERGED}" --input_type fastq --nproc "${THREADS}" \
    --bowtie2out "${OUT_DIR}/${SRA}_metaphlan_bowtie2.txt" \
    -o "${META}"
fi

# 3) HUMAnN run
echo "[HUMAnN] Functional profiling…"
humann \
  --input "${MERGED}" \
  --input-format fastq \
  --output "${OUT_DIR}" \
  --threads "${THREADS}" \
  --nucleotide-database "${HUMANN_CHOCOPHLAN}" \
  --protein-database    "${HUMANN_UNIREF}" \
  --search-mode uniref90 \
  --taxonomic-profile "${META}" \
  --remove-temp-output \
  --o-log "${OUT_DIR}/${SRA}_humann.log" \
  --log-level INFO

# 4) Post-processing (normalize, regroup, rename)
pushd "${OUT_DIR}" >/dev/null

humann_renorm_table --input ${SRA}_merged.clean_genefamilies.tsv \
                    --output ${SRA}_genefamilies_cpm.tsv --units cpm --update-snames

humann_renorm_table --input ${SRA}_merged.clean_pathabundance.tsv \
                    --output ${SRA}_pathabundance_rel.tsv --units relab --update-snames

humann_regroup_table --input ${SRA}_genefamilies_cpm.tsv --groups uniref90_ko  --output ${SRA}_KO_cpm.tsv
humann_regroup_table --input ${SRA}_genefamilies_cpm.tsv --groups uniref90_rxn --output ${SRA}_RXN_cpm.tsv

humann_rename_table  --input ${SRA}_KO_cpm.tsv  --names kegg       --output ${SRA}_KO_cpm_named.tsv
humann_rename_table  --input ${SRA}_RXN_cpm.tsv --names metacyc-rxn --output ${SRA}_RXN_cpm_named.tsv

popd >/dev/null

echo "[HUMAnN] Done → ${OUT_DIR}"
