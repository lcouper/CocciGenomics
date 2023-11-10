## Steps and scripts used to analyze Coccidiodides genomes 
Relevant code snippet for each shown below

#### 1. Index reference genome  
File name: bwamem_index
```
module load bwa-mem2/2.2.1
bwa-mem2 index CocciRefGenome.fna
```

#### 2. Align sequences to reference genome    
File name: bwamem
```
module load bwa-mem2/2.2.1
bwa-mem2 mem CocciRefGenome.fna SJV_1_1.fastq SJV_1_2.fastq > SJV_1_Aligned.sam
```

