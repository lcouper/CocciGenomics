#!/bin/bash
#SBATCH --job-name=CallVariants
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --output=CallVariants.out

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load bwa
module load  bcftools

bcftools mpileup --threads 12 -f /global/scratch/users/lcouper/SJV_Genomes/CocciRefGenome.fna -q 20 -Q 20 *.deduped.bam \
| bcftools call --threads 12 -mv --ploidy 1 -Oz -o VCFFILE

