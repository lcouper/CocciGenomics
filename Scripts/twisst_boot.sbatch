#!/bin/bash
#SBATCH --job-name=twisst_boot
#SBATCH --account=xxxxx
#SBATCH --partition=xxxxx
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100%25
#SBATCH --output=logs/twisst_boot_%A_%a.log

set -euo pipefail
module load anaconda3/2024.10-1-11.4 raxml/8.2.12
export PYTHONPATH="${PYTHONPATH:-}:$HOME/software/genomics_general"

BASE=$YOUR-PATH

SPEC="${1:-${SPEC:-}}"
if [[ -z "$SPEC" ]]; then
  echo "ERROR: no spec given." >&2
  echo "usage: sbatch $(basename "$0") <spec>" >&2
  echo "available specs in ${BASE}/groups:" >&2
  ls -1 "${BASE}"/groups/refs_*.tsv 2>/dev/null \
    | xargs -n1 basename 2>/dev/null | sed 's/^refs_//; s/\.tsv$//; s/^/  /' >&2
  exit 2
fi

#   python $HOME/software/genomics_general/VCF_processing/parseVCF.py \
#       -i ../cocci37.snps.vcf.gz --ploidy 1 --skipIndels -o cocci37.geno.gz
GENO=$BASE/cocci37.geno.gz
REALWIN=$BASE/cocci37.w100.data.tsv   # window structure from the unbootstrapped run
HELPERS=$BASE                      # boot_geno.py and twisst_summarise.py live beside the data

R=$SLURM_ARRAY_TASK_ID
REFS=$BASE/groups/refs_${SPEC}.tsv
WORK=$BASE/boot_${SPEC}/rep${R}    # own dir: RAxML scatters RAxML_* and .phy into cwd,
SUM=$BASE/boot_${SPEC}/summaries   # so concurrent tasks MUST NOT share a directory
mkdir -p "$WORK" "$SUM"
cd "$WORK"


### 0. checks on the spec ####################################################
[ -f "$REFS" ] || { echo "ERROR: $REFS not found"; exit 1; }
[ -f "$GENO" ] || { echo "ERROR: $GENO not found"; exit 1; }
awk 'NF!=2 {print "ERROR: line "NR" has "NF" fields: "$0; bad=1} END{exit bad+0}' "$REFS"
grep -qw Outgroup "$REFS" || { echo "ERROR: no Outgroup in $REFS"; exit 1; }

# reference group names come from the file, so nothing here is spec-specific
mapfile -t GRPS < <(awk '$2!="Outgroup"{print $2}' "$REFS" | sort -u)
GARGS=()
for g in "${GRPS[@]}"; do GARGS+=(-g "$g"); done
GARGS+=(-g Focal -g Outgroup)
echo "rep $R  spec $SPEC  groups: ${GRPS[*]}"

CLIN="Kern1 Kern2 Kern3 Kern4 Kern5 Kern6 Kern7 Kern8 Kern9 Kern10 Kern11 Kern12 Kern13 Kern14 Kern16 Kern17 Kern18 Kern19 Kern20 Kern21 Kern22 Kern23 Kern24 Kern25 Kern26 Kern27"
SOIL="22AC2 34B2 58B1 87A1 13B1 PS02PN14-1 157b2 L100 118a3 239a3b2"   # one genome per burrow
LEGACY=""                        # out-of-panel controls; set to "Coahuila_1 Guerrero_1 San_Diego" to restore
FOCAL="$CLIN $LEGACY $SOIL"      # trimming this is the main cost lever

# every focal and every reference must exist in the .geno, or twisst fails only AFTER the
# expensive tree step. Subshell disables pipefail so gzip's SIGPIPE from head does not abort.
HDR=$(set +o pipefail; gzip -cd "$GENO" | head -1)
for s in $FOCAL $(awk '{print $1}' "$REFS"); do
  grep -qw -- "$s" <<<"$HDR" || { echo "ERROR: '$s' not among the samples in $GENO"; exit 1; }
done

### 1. bootstrap the geno within each 100-SNP window ##########################
python "$HELPERS/boot_geno.py" -i "$GENO" -o boot.geno.gz -w 100 -s "$R"

### 2. rebuild every window tree from the resampled sites #####################
python "$HOME/software/genomics_general/phylo/raxml_sliding_windows.py" \
  -g boot.geno.gz -p boot.w100 \
  --windType sites -w 100 -M 50 --model GTRCAT \
  --raxml raxmlHPC-AVX -T "$SLURM_CPUS_PER_TASK" \
  --log boot.w100.raxlog.txt

# window structure must be identical across replicates, else they won't align.
# boot_geno.py keeps each block's original positions, so the count is deterministic.
NW_BOOT=$(( $(wc -l < boot.w100.data.tsv) - 1 ))
if [ -f "$REALWIN" ]; then
  NW_REAL=$(( $(wc -l < "$REALWIN") - 1 ))
  [ "$NW_BOOT" -eq "$NW_REAL" ] || { echo "ERROR: $NW_BOOT windows, expected $NW_REAL"; exit 1; }
fi
echo "rep $R: $NW_BOOT windows built"    # cross-check across reps afterwards, see header

### 3. weight every focal against the leave-one-out references ################
# References are clone-corrected (one genome per burrow), so dropping the focal is also
# dropping its whole clonal lineage -- no focal is ever matched against its own near-clone.
for f in $FOCAL; do
  awk -v f="$f" '$1!=f' "$REFS" > "groups_${f}.tsv"
  printf "%s\tFocal\n" "$f"     >> "groups_${f}.tsv"

  python "$HOME/software/twisst/twisst.py" \
    -t boot.w100.trees.gz \
    -w "w_${f}.csv.gz" \
    --groupsFile "groups_${f}.tsv" \
    "${GARGS[@]}" \
    --outgroup Outgroup --outputTopos "topos_${f}.txt" \
    2> "tw_${f}.log"
done

### 4. per-isolate scalars + per-window weights ###############################
python "$HELPERS/twisst_summarise.py" \
  -r "$R" -g "$REFS" \
  -s "$SUM/rep${R}.csv" -w "$SUM/rep${R}.windows.csv.gz" \
  w_*.csv.gz

### 5. clean up (keep reps 1-2 for inspection) ################################
cd "$BASE"
if [ "$R" -gt 2 ]; then rm -rf "$WORK"; fi
echo "rep $R done"
