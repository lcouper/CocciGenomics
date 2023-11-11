#!/bin/bash
#SBATCH --job-name=bwamem_align
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=SJV_Alignment.out

## commands to run:

cd /global/scratch/users/lcouper/SJV_Genomes/

module load bwa-mem2/2.2.1

for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
bwa-mem2 mem -t 12 CocciRefGenome.fna \
${base}_1.fastq ${base}_2.fastq > "${base}.aligned.sam"
done

