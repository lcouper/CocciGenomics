## Steps and scripts used to analyze Coccidiodides genomes 
Relevant code snippet for each shown below

#### 1. Index reference genome  
File name: bwamem_index
```
module load bwa-mem2/2.2.1
bwa-mem2 index CocciRefGenome.fna
```

#### 2. Align sequences to reference genome    
File name: bwamem_align
```
module load bwa-mem2/2.2.1
for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
bwa-mem2 mem -t 12 CocciRefGenome.fna \
${base}_1.fastq ${base}_2.fastq > "${base}.aligned.sam"
done
```

