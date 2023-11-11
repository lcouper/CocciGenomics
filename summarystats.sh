#!/bin/bash
#SBATCH --job-name=summarystats
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=summarystats.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load samtools/1.8

for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools flagstat ${base}.bam > "${base}.bam.stats.txt"
done
