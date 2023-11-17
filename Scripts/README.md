## Steps and scripts used to analyze Coccidiodides genomes 
Relevant code snippet for each shown below

#### 1. Index reference genome  
Script name: bwamem_index
```
bwa-mem2 index CocciRefGenome.fna
```

#### 2. Align sequences to reference genome    
Script name: bwamem_align
```
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
for infile in Aligned/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b Aligned/${base}.aligned.sam > Aligned/"${base}.aligned.bam"
done
```

#### 4. Sort bam file by coordinates using samtools
Script name: sortbam.sh
```
for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools sort -@ 12 -o "${base}.sorted.bam" ${base}.bam
done
```

### 5. Obtain summary stats about bam file
Script name: summarystats.sh
```
for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools flagstat ${base}.bam > "${base}.bam.stats.txt"
done
```

### 6. Mark and remove duplicates using picard
Note: picard.jar was downloaded from the Broad institute here: https://github.com/broadinstitute/picard/releases/tag/3.1.1
Script name: MarkDups.sh
```
for infile in *.aligned.sorted.bam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sorted.bam)
java -jar picard.jar MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I ${base}.aligned.sorted.bam \
-O "${base}.deduped.bam" \
-M "${base}.dup_metrics.txt"
done
```
