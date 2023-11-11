## Steps and scripts used to analyze Coccidiodides genomes 
Relevant code snippet for each shown below

#### 1. Index reference genome  
Script name: bwamem_index
```
module load bwa-mem2/2.2.1
bwa-mem2 index CocciRefGenome.fna
```

#### 2. Align sequences to reference genome    
Script name: bwamem_align
```
module load bwa-mem2/2.2.1
for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
bwa-mem2 mem -t 12 CocciRefGenome.fna \
${base}_1.fastq ${base}_2.fastq > "${base}.aligned.sam"
done
```

#### 3. Compress .sam to .bam using samtools
Script name: sam2bam.sh
```
module load samtools/1.8

for infile in Aligned/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b Aligned/${base}.aligned.sam > Aligned/"${base}.aligned.bam"
done
```
