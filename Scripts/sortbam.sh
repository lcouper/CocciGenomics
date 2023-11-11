#!/bin/bash
#SBATCH --job-name=sortbam
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=sortbam.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned

module load samtools/1.8

for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools sort -@ 12 -o "${base}.sorted.bam" ${base}.bam
done

