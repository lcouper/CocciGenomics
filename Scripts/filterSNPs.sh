#!/bin/bash
#SBATCH --job-name=filterSNPS
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=filterSNPs.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams

module load bamtools/2.5.1

bcftools view -O z -o Filtered_VCFFILE.vcf.gz -e 'QUAL<=40' VCFFILE.vcf
