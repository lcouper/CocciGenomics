## Steps and scripts used to process WGS data on environmetnal Coccidiodides samples ### 

Relevant code snippet for each shown below

### 1. Obtained raw reads from Berkeley QB3.  
The fastq.gz files (1 forward, 1 reverse) are stored here on the Remais Group Shared Drive: SPORE/WGS/Sequence data (All)/ 

### 2. Filter poor quality reads and trim poor quality bases 

Software used: Trimmomatic V 0.39 (Bolger et al. 2014)
Script: trim.sbatch
Code snippet for single sample:
```
module load bio/trimmomatic/0.39-gcc-11.4.0
trimmomatic PE PS02PN14-1_S1_L007_R1_001.fastq.gz PS02PN14-1_S1_L007_R2_001.fastq.gz \
PS02PN14-1_S1_L007_R1_001.trim.fastq.gz PS02PN14-1_S1_L007_R1_001.untrim.fastq.gz \
PS02PN14-1_S1_L007_R2_001.trim.fastq.gz PS02PN14-1_S1_L007_R2_001.untrim.fastq.gz\
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```
### 3. 

Software used: bio/fastqc/0.12.1-gcc-11.4.0
Script: fastqc.sbatch
Code snippet for single sample:
```
module load bio/fastqc/0.12.1-gcc-11.4.0
fastqc trimmed_fastqc/*.fastq*
```

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



