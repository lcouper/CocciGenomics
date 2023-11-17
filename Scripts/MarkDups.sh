#!/bin/bash
#SBATCH --job-name=MarkDups
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=08:00:00
#SBATCH --output=MarkDups.out

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load java
module load bwa

for infile in *.aligned.sorted.bam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sorted.bam)
java -jar picard.jar MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I ${base}.aligned.sorted.bam \
-O "${base}.deduped.bam" \
-M "${base}.dup_metrics.txt"
done
