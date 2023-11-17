#!/bin/bash
#SBATCH --job-name=IndexBam
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=08:00:00
#SBATCH --output=IndexBam.out

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load samtools

for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools index -b ${base}.bam
done
