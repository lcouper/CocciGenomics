## Steps and scripts used to process WGS data on environmetnal _Coccidiodides_ samples ### 

Relevant code snippet for each shown below

### 1. Obtained raw reads from Berkeley QB3.  

The fastq.gz files (1 forward, 1 reverse) are stored here on the Remais Group Shared Drive: SPORE/WGS/Sequence data (All)/ 
and Berkeley's HPC BRC at: /global/scratch/users/lcouper/SoilCocciSeqs

### 2. Filter poor quality reads and trim poor quality bases 

Software used: Trimmomatic V 0.39 (Bolger et al. 2014)   
Script: trim.sh
Code snippet for single sample:    
```
module load bio/trimmomatic/0.39-gcc-11.4.0
trimmomatic PE PS02PN14-1_S1_L007_R1_001.fastq.gz PS02PN14-1_S1_L007_R2_001.fastq.gz \
PS02PN14-1_S1_L007_R1_001.trim.fastq.gz PS02PN14-1_S1_L007_R1_001.untrim.fastq.gz \
PS02PN14-1_S1_L007_R2_001.trim.fastq.gz PS02PN14-1_S1_L007_R2_001.untrim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```

### 3. Perform quality check on samples using fastqc

Software used: bio/fastqc/0.12.1-gcc-11.4.0   
Script: fastqc.sh    
Code snippet for single sample:      
```
module load bio/fastqc/0.12.1-gcc-11.4.0
fastqc trimmed_fastqc/*.fastq.gz
```

### 4. Index reference genome  

Note: Using reference genome for Coccidioides immitis RS (GCA_000149335.2)   
Downloaded here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/   
Saved/uploaded as: CocciRef_GCA_000149335.2.fna  

Software used: bio/bwa-mem2/2.2.1  
Script name: bwamem_index.sh    
Code snippet:   

```
module load bio/bwa-mem2/2.2.1
bwa-mem2 index CocciRef_GCA_000149335.2.fna
```

### 5. Align sequences to reference genome    

Software used: bio/bwa-mem2/2.2.1   
Script name: alignreads.sh    
Code snippet:   

```
#First unzip trimmed fastq files if not done already
gunzip trimmed_fastq/*.gz

#Then align to cocci immitis RS ref genome

for infile in trimmed_fastq/*_R1_001.trim.fastq
do
base=$(basename ${infile} _R1_001.trim.fastq)
bwa-mem2 mem -t 12 RefGenome/CocciRef_GCA_000149335.2.fna \
trimmed_fastq/${base}_R1_001.trim.fastq trimmed_fastq/${base}_R2_001.trim.fastq > results/sam/"${base}.aligned.sam"
done
```

### 6. Compress sam to bam, sort bam files, and extract mapping stats

Software used: bio/samtools/1.17-gcc-11.4.0    
Script name: sam2bam.sh      
Code snippet:     

```
for infile in results/sam/*.aligned.sam
do
echo "working with file $infile"
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b results/sam/${base}.aligned.sam > results/bam/"${base}.aligned.bam"
samtools sort -@ 12 results/bam/${base}.aligned.bam -o "${base}.sorted.bam"
samtools flagstat results/bam/${base}.aligned.bam > results/bam/"${base}.bam.stats.txt"
done
```

### 7. Mark and remove duplicates 

Software used: bio/picard/3.0.0-gcc-11.4.0       
Script name: markdups.sh    
Code snippet:

```
for infile in results/sortedbams/*.sorted.bam
do
base=$(basename ${infile} .sorted.bam)
picard MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I results/sortedbams/${base}.sorted.bam \
-O results/dedupedbams/"${base}.deduped.bam" \
-M results/sortedbams/"${base}.dup_metrics.txt"
done
```

### 8. Index de-duplicated bam files 

Software used: bio/samtools/1.17-gcc-11.4.0    
Script name: indexbam.sh
Code snippet:

```
for infile in results/dedupedbams/*.deduped.bam
do
base=$(basename ${infile} .deduped.bam)
samtools index -b results/dedupedbams/${base}.deduped.bam
done
```

### 9. Compute alignment statistics

Note: calculates statistics including total reads, mapped reads, % failed QC, % duplicates, % paired-end reads, % singletons

Software used: bio/bamtools/2.5.2-gcc-11.4.0
Script: alignstats.sh  
Code snippet:

```
for infile in results/dedupedbams/*.deduped.bam
do
base=$(basename ${infile} .deduped.bam)
bamtools stats -in results/dededupedbams/${base}.deduped.bam > results/dedupedbams/"${base}.AlignStats.txt"
done
```

#### Look into repeat masking! 


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

### 14. Filter SNVs using vcftools and remove multi-allelic sites 
Script name: vcf_filter.sh

* Note vcftools installed by BRC into the module farm. Load and run using the following:
```
module unload gcc/6.3.0
module load gcc/11.3.0
export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles
module load vcftools/0.1.16
module load bcftools/1.6

vcftools --vcf Filtered_Sorted_VCFFILE.vcf --maf 0.05 --minQ 30 --max-missing 0.75 --minDP 10 --recode --recode-INFO-all --out VCF_AllVariants.vcf
# 202,550 out of a possible 260,032 Sites retained

bcftools view -m2 -M2 -v snps VCF_AllVariants.vcf > VCF_Biallelic.vcf

# 194668 SNPs retained
```
To identify number of SNPs in vcf file:
grep -v "^#" VCF_Biallelic.vcf|wc -l

* Note: 
uploaded VCFfile to SCG (to ThermalSelectionExpSeqFiles > results > bam > deduped_bams > filtered_VCF. Name “Filtered_Sorted_VCFFILE_SJV_Genomes.vcf”

### 15. Generate genotype matrix
Done using vcftools. Outputs 3 files: ‘.012’ contains the genotypes of each individual on a separate line (with 0, 1, 2 denoting the number of non-reference alleles at the site), ‘.ind’ lists the individuals included in the main file, ‘.pos’ details the site location included in the main file. 
```
vcftools --012 --vcf Filtered_Sorted_VCFFILE_SJV_Genomes.vcf --out SJV_genotype_matr
```

## Additional downstream steps of interest:

### 1. calculate # of SNPs differing between each possible pair using plink.
Note I run this on SCG, using the pre-filtered VCFFILE. Followed guidance here: https://www.biostars.org/p/351404/
Note: Prior work found 1000s of SNPs separate strains from the same region  (mentioned here: https://academic.oup.com/cid/article/60/1/e1/2895394)

```
module load plink/1.90
plink --vcf VCF_Biallelic.vcf --allow-extra-chr --genome full --out plink.genome.SJV
# in resulting --out plink.genome.SJV.genome file, if you sum the IBS0 and IBS1 columns (in R), you get the pairwise SNP differences
```

### 2. Generate phylogenetic tree
Note this is done in R, using the genotype matrix created above (i.e., SJV_genotype_matr) and the ape package
```
SJV_gm <- fread("GenotypeMatrix/SJV_genotype_matr.csv", header = TRUE)
stree = nj(dist.gene(datasub)) # this specifies the 'neighbor-joining' method for tree-construction (appears to be most common)


### 3. Plot PCA
Note this is done in R, using the genotype matrix created above (i.e., SJV_genotype_matr) along with metadata labels



