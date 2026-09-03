#!/bin/bash
#SBATCH --job-name=mat_diagnostic_assembly
#SBATCH --account=xxxxxx
#SBATCH --partition=yyyyyy
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=mat_diagnostic_assembly.out

set -euo pipefail
shopt -s nullglob

module load bio/bwa-mem2/2.2.1
module load bio/samtools/1.17-gcc-11.4.0
module load bio/blast-plus/2.13.0-gcc-11.4.0
module load python/3.10.12-gcc-11.4.0
module load spades/4.1.0

echo "Python: $(which python)"
python --version
echo "SPAdes: $(which spades.py)"
spades.py --version

FASTQ_DIR=$YOUR-FASTQ-DIR
MAT_DIR=$YOUR-MAT-DIR
OUT_DIR=${MAT_DIR}/MAT_diagnostic_assemblies_from_fastq

MAT11_FULL=${MAT_DIR}/EF472259.1.fna
MAT12_FULL=${MAT_DIR}/EF472258.1.fna

MAT11_DIAG=${MAT_DIR}/MAT1-1_diag_EF472259_10824_12046.fna
MAT12_DIAG=${MAT_DIR}/MAT1-2_diag_EF472258_1_3289.fna
MAT_DIAG_REF=${MAT_DIR}/MAT_diagnostic_refs.fna
MAT_DIAG_DB=${MAT_DIR}/MAT_diagnostic_refs_DB

DETAIL_SUMMARY=${OUT_DIR}/MAT_diagnostic_blastn_summary.tsv
CALL_SUMMARY=${OUT_DIR}/MAT_diagnostic_calls.tsv
READCOUNT_SUMMARY=${OUT_DIR}/MAT_diagnostic_read_counts.tsv

THREADS=8
MQ_MIN=20
MIN_BAITED_READS=10
MIN_PID=90
MIN_QCOV=80
EXCL_FLAGS=$((4|256|2048))

mkdir -p "$OUT_DIR"

for f in "$MAT11_FULL" "$MAT12_FULL"; do
  [ -s "$f" ] || { echo "Missing required file: $f" >&2; exit 1; }
done

echo "Creating diagnostic MAT reference files in: $MAT_DIR"

samtools faidx "$MAT11_FULL"
samtools faidx "$MAT12_FULL"

samtools faidx "$MAT11_FULL" EF472259.1:10824-12046 | sed '1s/.*/>MAT1-1_diag_EF472259_10824_12046/' > "$MAT11_DIAG"
samtools faidx "$MAT12_FULL" EF472258.1:1-3289 | sed '1s/.*/>MAT1-2_diag_EF472258_1_3289/' > "$MAT12_DIAG"

cat "$MAT11_DIAG" "$MAT12_DIAG" > "$MAT_DIAG_REF"

echo "Diagnostic reference lengths:"
awk '/^>/ {if (seqlen) print seqlen; print; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$MAT_DIAG_REF"

echo "Checking for cross-hits between diagnostic references:"
makeblastdb -in "$MAT_DIAG_REF" -dbtype nucl -out "$MAT_DIAG_DB" >/dev/null
blastn -query "$MAT_DIAG_REF" -db "$MAT_DIAG_DB" -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" | awk '$1 != $2' | sort -k12,12nr -k4,4nr || true

bwa-mem2 index "$MAT_DIAG_REF" >/dev/null 2>&1 || true

echo -e "Sample\tQuery\tqseqid\tsseqid\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tN_contigs\tBaited_reads\tReads_MAT1_1\tReads_MAT1_2" > "$DETAIL_SUMMARY"
echo -e "Sample\tCall\tReason\tBaited_reads\tReads_MAT1_1\tReads_MAT1_2\tBest_MAT1_1_pid\tBest_MAT1_1_length\tBest_MAT1_1_qcov\tBest_MAT1_2_pid\tBest_MAT1_2_length\tBest_MAT1_2_qcov\tN_contigs" > "$CALL_SUMMARY"
echo -e "Sample\tReference\tMapped_reads" > "$READCOUNT_SUMMARY"

for R1 in ${FASTQ_DIR}/*_1.len75.trim.fastq; do
  SAMPLE=$(basename "$R1" _1.len75.trim.fastq)
  R2=${FASTQ_DIR}/${SAMPLE}_2.len75.trim.fastq

  if [ ! -s "$R2" ]; then echo "Missing R2 for ${SAMPLE}: ${R2}" >&2; continue; fi

  SAMPLE_DIR=${OUT_DIR}/${SAMPLE}
  mkdir -p "$SAMPLE_DIR"
  cd "$SAMPLE_DIR"

  echo "=== ${SAMPLE} ==="

  bwa-mem2 mem -t "$THREADS" "$MAT_DIAG_REF" "$R1" "$R2" | samtools view -b -F "$EXCL_FLAGS" -q "$MQ_MIN" -o ${SAMPLE}.diagnostic_MAT.bam

  BAITED_READS=$(samtools view -c ${SAMPLE}.diagnostic_MAT.bam || echo 0)
  READS_MAT11=$(samtools view ${SAMPLE}.diagnostic_MAT.bam | awk '$3=="MAT1-1_diag_EF472259_10824_12046" {n++} END {print n+0}')
  READS_MAT12=$(samtools view ${SAMPLE}.diagnostic_MAT.bam | awk '$3=="MAT1-2_diag_EF472258_1_3289" {n++} END {print n+0}')

  echo "  baited diagnostic MAT reads: ${BAITED_READS}"
  echo "  MAT1-1 diagnostic reads: ${READS_MAT11}"
  echo "  MAT1-2 diagnostic reads: ${READS_MAT12}"

  echo -e "${SAMPLE}\tMAT1-1_diag_EF472259_10824_12046\t${READS_MAT11}" >> "$READCOUNT_SUMMARY"
  echo -e "${SAMPLE}\tMAT1-2_diag_EF472258_1_3289\t${READS_MAT12}" >> "$READCOUNT_SUMMARY"

  if [ "$BAITED_READS" -lt "$MIN_BAITED_READS" ]; then
    echo -e "${SAMPLE}\tLOW_READS\tBaited reads below ${MIN_BAITED_READS}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}\t0\t0\t0\t0\t0\t0\t0" >> "$CALL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-1\tLOW_READS\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-2\tLOW_READS\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo "  → skipped assembly; too few baited reads"
    continue
  fi

  samtools sort -n -@ "$THREADS" -o ${SAMPLE}.diagnostic_MAT.namesort.bam ${SAMPLE}.diagnostic_MAT.bam
  samtools fastq -@ "$THREADS" -1 ${SAMPLE}.diagnostic_MAT.R1.fq -2 ${SAMPLE}.diagnostic_MAT.R2.fq -s ${SAMPLE}.diagnostic_MAT.single.fq -0 /dev/null -N ${SAMPLE}.diagnostic_MAT.namesort.bam

  if ! spades.py -1 ${SAMPLE}.diagnostic_MAT.R1.fq -2 ${SAMPLE}.diagnostic_MAT.R2.fq -s ${SAMPLE}.diagnostic_MAT.single.fq -o spades_diagnostic_mat -t "$THREADS" --only-assembler --cov-cutoff off -k 21,33 > spades_diagnostic.stdout 2> spades_diagnostic.stderr; then
    echo -e "${SAMPLE}\tSPADES_FAILED\tSPAdes failed\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}\t0\t0\t0\t0\t0\t0\t0" >> "$CALL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-1\tSPADES_FAILED\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-2\tSPADES_FAILED\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo "  → SPAdes failed; continuing to next sample"
    continue
  fi

  CONTIGS=spades_diagnostic_mat/contigs.fasta

  if [ ! -s "$CONTIGS" ]; then
    echo -e "${SAMPLE}\tNO_CONTIGS\tNo SPAdes contigs produced\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}\t0\t0\t0\t0\t0\t0\t0" >> "$CALL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-1\tNO_CONTIGS\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo -e "${SAMPLE}\tMAT1-2\tNO_CONTIGS\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t0\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
    echo "  → no contigs produced"
    continue
  fi

  N_CONTIGS=$(grep -c "^>" "$CONTIGS")

  blastn -query "$CONTIGS" -db "$MAT_DIAG_DB" -evalue 1e-20 -perc_identity 80 -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" > ${SAMPLE}.diagnostic_contig_blast.tsv || true

  HIT11=$(awk '$2=="MAT1-1_diag_EF472259_10824_12046"' ${SAMPLE}.diagnostic_contig_blast.tsv | sort -k12,12nr -k4,4nr | head -n1)
  HIT12=$(awk '$2=="MAT1-2_diag_EF472258_1_3289"' ${SAMPLE}.diagnostic_contig_blast.tsv | sort -k12,12nr -k4,4nr | head -n1)

  if [ -n "$HIT11" ]; then echo -e "${SAMPLE}\tMAT1-1\t${HIT11}\t${N_CONTIGS}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
  else echo -e "${SAMPLE}\tMAT1-1\tNO_HIT\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t${N_CONTIGS}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"; fi

  if [ -n "$HIT12" ]; then echo -e "${SAMPLE}\tMAT1-2\t${HIT12}\t${N_CONTIGS}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"
  else echo -e "${SAMPLE}\tMAT1-2\tNO_HIT\tNA\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0\t${N_CONTIGS}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}" >> "$DETAIL_SUMMARY"; fi

  MAT11_STATS=$(awk -v minpid="$MIN_PID" -v minqcov="$MIN_QCOV" '$2=="MAT1-1_diag_EF472259_10824_12046" {qcov=100*$4/$6; if ($3>=minpid && qcov>=minqcov && $12>best) {best=$12; pid=$3; len=$4; cov=qcov}} END {printf "%s\t%s\t%.2f\t%s", pid+0, len+0, cov+0, (pid>=minpid && cov>=minqcov ? 1 : 0)}' ${SAMPLE}.diagnostic_contig_blast.tsv)
  MAT12_STATS=$(awk -v minpid="$MIN_PID" -v minqcov="$MIN_QCOV" '$2=="MAT1-2_diag_EF472258_1_3289" {qcov=100*$4/$6; if ($3>=minpid && qcov>=minqcov && $12>best) {best=$12; pid=$3; len=$4; cov=qcov}} END {printf "%s\t%s\t%.2f\t%s", pid+0, len+0, cov+0, (pid>=minpid && cov>=minqcov ? 1 : 0)}' ${SAMPLE}.diagnostic_contig_blast.tsv)

  MAT11_PID=$(echo "$MAT11_STATS" | cut -f1)
  MAT11_LEN=$(echo "$MAT11_STATS" | cut -f2)
  MAT11_QCOV=$(echo "$MAT11_STATS" | cut -f3)
  MAT11_PASS=$(echo "$MAT11_STATS" | cut -f4)
  MAT12_PID=$(echo "$MAT12_STATS" | cut -f1)
  MAT12_LEN=$(echo "$MAT12_STATS" | cut -f2)
  MAT12_QCOV=$(echo "$MAT12_STATS" | cut -f3)
  MAT12_PASS=$(echo "$MAT12_STATS" | cut -f4)

  if [ "$MAT11_PASS" -eq 1 ] && [ "$MAT12_PASS" -eq 0 ]; then CALL="MAT1-1"; REASON="Diagnostic MAT1-1 contig hit passes thresholds"
  elif [ "$MAT11_PASS" -eq 0 ] && [ "$MAT12_PASS" -eq 1 ]; then CALL="MAT1-2"; REASON="Diagnostic MAT1-2 contig hit passes thresholds"
  elif [ "$MAT11_PASS" -eq 1 ] && [ "$MAT12_PASS" -eq 1 ]; then CALL="AMBIGUOUS_BOTH"; REASON="Both diagnostic contig hits pass thresholds; inspect manually"
  else CALL="NO_DIAGNOSTIC_CONTIG_HIT"; REASON="No diagnostic contig hit passes thresholds"; fi

  echo -e "${SAMPLE}\t${CALL}\t${REASON}\t${BAITED_READS}\t${READS_MAT11}\t${READS_MAT12}\t${MAT11_PID}\t${MAT11_LEN}\t${MAT11_QCOV}\t${MAT12_PID}\t${MAT12_LEN}\t${MAT12_QCOV}\t${N_CONTIGS}" >> "$CALL_SUMMARY"

  awk -v s="$SAMPLE" '/^>/ {sub(/^>/, ">" s "|"); print; next} {print}' "$CONTIGS" > ${SAMPLE}.diagnostic_MAT_contigs.fasta

  echo "  → call: ${CALL}; contigs: ${N_CONTIGS}"
done

echo
echo "Done."
echo "Call summary: $CALL_SUMMARY"
echo "Detailed BLAST summary: $DETAIL_SUMMARY"
echo "Read-count summary: $READCOUNT_SUMMARY"
