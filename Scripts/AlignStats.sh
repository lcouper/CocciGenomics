#!/bin/bash
#SBATCH --job-name=AlignStats
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=12
#SBATCH --time=08:00:00
#SBATCH --output=AlignStats.out

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load bamtools

for infile in *.deduped.bam
do
echo "working with file $infile"
base=$(basename ${infile} .deduped.bam)
bamtools stats -in ${base}.deduped.bam > "${base}.AlignStats.txt"
done
