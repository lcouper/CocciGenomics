## Steps and scripts used to analyze Coccidiodides genomes 
Relevant code snippet for each shown below

### 1. Index reference genome  
Script name: bwamem_index
```
bwa-mem2 index CocciRefGenome.fna
```

### 2. Align sequences to reference genome    
Script name: bwamem_align
```
for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
bwa-mem2 mem -t 12 CocciRefGenome.fna \
${base}_1.fastq ${base}_2.fastq > "${base}.aligned.sam"
done
```

### 3. Compress .sam to .bam using samtools
Script name: sam2bam.sh
```
for infile in Aligned/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b Aligned/${base}.aligned.sam > Aligned/"${base}.aligned.bam"
done
```

### 4. Sort bam file by coordinates using samtools
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

### 6. Mark and remove duplicates 
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

### 7. Index de-duplicated bam files 
Script name: IndexBam.sh
```
for infile in *.bam
do
echo "working with file $infile"
base=$(basename ${infile} .bam)
samtools index -b ${base}.bam
done
```

### 8. Compute alignment statistics
Note: calculates statistics including total reads, mapped reads, % failed QC, % duplicates, % paired-end reads, % singletons
Script: AlignStats.sh   
```
for infile in *.deduped.bam
do
echo "working with file $infile"
base=$(basename ${infile} .deduped.bam)
bamtools stats -in ${base}.deduped.bam > "${base}.AlignStats.txt"
done
```

### 9. Detect single nucleotide variants
Script: CallVariants.sh
```
bcftools mpileup --threads 12 -f /global/scratch/users/lcouper/SJV_Genomes/CocciRefGenome.fna -q 20 -Q 20 *.deduped.bam \
| bcftools call --threads 12 -mv --ploidy 1 -Oz -o VCFFILE
```

### 10. Filter low quality variants
Script: FilterSNP.sh
```
bcftools view -O z -o Filtered_VCFFILE.vcf.gz -e 'QUAL<=40' VCFFILE.vcf
```

### 11. Index filtered vcffiles
* Note: required installing tabix. Followed instructions here : https://github.com/trinityrnaseq/Griffithlab_rnaseq_tutorial_wiki/blob/master/AWS-Setup.md to install to SJV_Genomes/Aligned/SortedBams
* Set export path as: export PATH=$PATH:/global/scratch/users/lcouper/SJV_Genomes/Aligned/SortedBams/tabix-0.2.6
```
tabix *.vcf
```

### 12. Create dictionary for reference genome
* Note picard.jar uploaded to working directory
```
module load java
java -jar picard.jar CreateSequenceDictionary -R ../../CocciRefGenome.fna -O ../../CocciRef.dict
```

### 13. Sort vcf according to refeference dictionary 
```
module load bwa
module load java
java -jar picard.jar SortVcf \
-I Filtered_VCFFILE.vcf.gz \
-O Filtered_Sorted_VCFFILE.vcf.gz \
-SD /global/scratch/users/lcouper/SJV_Genomes/CocciRef.dict
```

### 14. Filter SNVs using vcftools
* will require installing vcftools. Have requested this from BCR help


## Additional downstream steps of interest:
- calculate # of SNPs differing between each possible pair (can maybe be done with program 'plink' (Available as a module on savio) according to : https://www.biostars.org/p/351404/ 
- generate phylogenetic tree: have seen 'MrBayes' as gold standard for this, but have never used it before
- calculate overall nucleotide diversity
- plot PCA

Downloaded MrBayes here: https://github.com/NBISweden/MrBayes/tree/v3.2.7a
followed instruction manual here: https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf
may require using: https://github.com/edgardomortiz/vcf2phylip to get VCFFILE into Nexus format for use with MrBayes

