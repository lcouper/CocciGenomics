#!/bin/bash
#SBATCH --job-name=sam2bam
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=sam2bam.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/

module load samtools/1.8

for infile in Aligned/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b Aligned/${base}.aligned.sam > Aligned/"${base}.aligned.bam"
done

