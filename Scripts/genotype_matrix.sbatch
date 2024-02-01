#!/bin/bash
#SBATCH --job-name=genotypeMatrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=100GB
#SBATCH --partition=batch
#SBATCH --account=emordeca
#SBATCH --time=14:00:00
#SBATCH --output=genotype_matrix_SJV.out

module load bwa
module load bcftools
module load vcftools

cd /labs/emordeca/ThermalSelectionExpSeqFiles/results/bam/deduped_bams/filtered_VCF

vcftools --012 --vcf Filtered_Sorted_VCFFILE_SJV_Genomes.vcf --out SJV_genotype_matrix

