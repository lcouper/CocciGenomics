#!/bin/bash
#SBATCH --job-name=test
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=00:10:00
#SBATCH --output=SJV_1_Align.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/

module load bwa-mem2/2.2.1

bwa-mem2 index CocciRefGenome.fna
