## Steps and scripts used to process _Coccidiodides_ genomes from environmental samples or SRA ### 

Relevant code snippet for each shown below

### 1. Obtained raw reads from Berkeley QB3 or SRA.  

The fastq.gz files (1 forward, 1 reverse) are stored here on the Remais Group Shared Drive: SPORE/WGS/Sequence data (All)/ 
and Berkeley's HPC BRC at: /global/scratch/users/lcouper/SoilCocciSeqs

Prior cocci sequences were downloaded from NCBI using the SRA toolkit.   
Additional notes [here](https://docs.google.com/document/d/1gkM7m6TjQAOO1pwxe4X2DrIuMPuA3uGd6UalImpb-h4/edit?tab=t.0).   
Tracker for sequences downloaded and metadata [here](https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit?gid=1963297784#gid=1963297784). 

### 2. Filter poor quality reads and trim poor quality bases 

Note that Illumina adapters [available and downlaoded from here](https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE.fa)    
Software used: Trimmomatic V 0.39 (Bolger et al. 2014)      
Job Script: trim.sh and trim.sra.sh   
Relevant code snippet:       
```
module load bio/trimmomatic/0.39-gcc-11.4.0
trimmomatic PE PS02PN14-1_S1_L007_R1_001.fastq.gz PS02PN14-1_S1_L007_R2_001.fastq.gz \
PS02PN14-1_S1_L007_R1_001.trim.fastq.gz PS02PN14-1_S1_L007_R1_001.untrim.fastq.gz \
PS02PN14-1_S1_L007_R2_001.trim.fastq.gz PS02PN14-1_S1_L007_R2_001.untrim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```

### 3. Perform quality check on samples using fastqc

Software used: bio/fastqc/0.12.1-gcc-11.4.0   
Script: fastqc.sh, fastqc.sra.sh        
Relevant code snippet:      
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
Script name: alignreads.sh, alignreads.sra.sh    
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
Script name: sam2bam.sh, sam2bam.sra.sh          
Code snippet:     

```
for infile in results/sam/*.aligned.sam
do
base=$(basename ${infile} .aligned.sam)
samtools view -@ 12 -b results/sam/${base}.aligned.sam > results/bam/"${base}.aligned.bam"
samtools sort -@ 12 results/bam/${base}.aligned.bam -o "${base}.sorted.bam"
samtools flagstat results/bam/${base}.aligned.bam > results/bam/"${base}.bam.stats.txt"
done
```

### 7. Mark and remove duplicates 

Software used: bio/picard/3.0.0-gcc-11.4.0       
Script name: markdups.sh, markdups.sra.sh,  
Code snippet:

```
for infile in results/sortedbams/*.sorted.bam
do
base=$(basename ${infile} .sorted.bam)
picard MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I results/sortedbams/${base}.sorted.bam \
-O results/dedupedbams/"${base}.deduped.bam" \
-M results/dedupedbams/"${base}.dup_metrics.txt"
done
```

## Option 1: Use GATK for variant calling 

### 8. Add read groups 

*Followed guidance here: https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information
and issue was diagnosed here: https://gatk.broadinstitute.org/hc/en-us/community/posts/4412745467931-HaplotypeCaller-does-not-work

See this spreadsheet for what read group parameters were added:
https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit?gid=1963297784#gid=1963297784


Software used: bio/picard/3.0.0-gcc-11.4.0, java  
Script name: addrg.sbatch
Code snippet:

```
module load java
module load bio/picard/3.0.0-gcc-11.4.0 

picard AddOrReplaceReadGroups \
I=results/dedupedbams/PS02PN14-2_S2_L007.deduped.bam \
O=results/bamswithrg/PS02PN14-2_S2_L007.rg.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=21
```

### 9. Index bam files with read group added

Software used: bio/samtools/1.17-gcc-11.4.0    
Code snippet:

```
samtools index -b results/bamswithrg/${base}.rg.bam
done
```

### 10. Call variants using GATK HaplotypeCaller 

##### 10a. Call variants #####
*Note that I downloaded gatk from [here](https://github.com/broadinstitute/gatk/releases) and then uploaded the jar file to savio to the SoilCocciSeqs directory
Guidance on these steps found here: https://www.biostars.org/p/405702/   
Software used: java, gatk 4.5.0.0
Script name: haplo.sh
Code snippet:

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" HaplotypeCaller \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
-ploidy 1 \
-ERC BP_RESOLUTION \
-I results/bamswithrg/B0727_Argentina.rg.bam \    # Note this code is run on each sample individually (could also do in loop, but does take a while to run for each)
--output-mode EMIT_ALL_CONFIDENT_SITES \
-O results/haplocalled/B0727_Argentina.g.vcf.gz
```

### 10b. Combine GVCF files 

*Note, combined all the above files into a single directory 'AllGenomesHaploCalled'. Then, created a list of files in this directory using*
```
ls AllGenomesHaploCalled/*.vcf.gz > gvcfs.list
```

Software used: java, gatk 4.5.0.0
Script name: combinegvcfs.sh

```
module load java

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" CombineGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
--variant gvcfs.list \
-O combined.g.vcf.gz

```

### 10c. Joint-genotyping on combined GVCF files 

Software used: java, gatk 4.5.0.0
Script name: genotypegvcfs.sh

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" GenotypeGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
-ploidy 1 \
# --include-non-variant-sitess \    <- note this is an optional argument for retaining all SNPs, which I did not want/use
-V combined.g.vcf.gz \
-O final.withoutNonSNPs.vcf.gz
```

Next, unzip final.withoutNonSNPs.vcf.gz file and identify number of variant sites:

```
gunzip final.withoutNonSNPs.vcf.gz
grep -v "^#" AllGenomesHaploCalled/final.withoutNonSNPs.vcf | wc -l   # 748,817
```


### 10d. Filter variants 

Software used: java, gatk 4.5.0.0    
Script: filtervcfs.sbatch    
Code snippet:    

```
module load java

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" VariantFiltration \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
--variant final.withoutNonSNPs.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AC > 2 || DP < 10" \
--filter-name "AllFilters" \
-O final.filtered.withoutNonSNPs.vcf

```



### 103. Select variants 

Software used: java, gatk 4.5.0.0    
Script: selectsnps.sbatch    
Code snippet:    

```
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SelectVariants \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
--variant final.filtered.withoutNonSNPs.vcf \
--select-type SNP \
-O final.SNPs.vcf # Note that 683,199 SNPs remain 
```




## Option 2: Use bcftools for variant calling (this one giving issues 

### 10. Index de-duplicated bam files 

Software used: bio/samtools/1.17-gcc-11.4.0    
Script name: indexbam_dedupedbams.sh, indexbam_dedupedbams.sra.sh    
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
Script: alignstats.sh, alignstats.sra.sh     
Code snippet:  

```
for infile in results/dedupedbams/*.deduped.bam
do
base=$(basename ${infile} .deduped.bam)
bamtools stats -in results/dededupedbams/${base}.deduped.bam > results/dedupedbams/"${base}.AlignStats.txt"
done
```

### 9.1 Compute depth at each position of sample 

Software used: bio/samtools/1.17-gcc-11.4.0    
Script name: depth.sh, depth.sra.sh     
Code snippet:

```
for infile in *.bam
do
base=$(basename ${infile} .bam)
samtools depth -a ${base}.bam > "${base}.depth.txt"
done
```

### 10. Detect single nucleotide variants
Script: callvariants.sh
```
bcftools mpileup --threads 12 -f RefGenome/CocciRef_GCA_000149335.2.fna -q 20 -Q 20 results/dedupedbams/*.deduped.bam \
| bcftools call --threads 12 -mv --ploidy 1 -Oz -o VCFFILE
```

### 10. Filter low quality variants

Software used: bio/bcftools/1.16-gcc-11.4.0
Script: filtervariants.sh

```
bcftools view -O z -o results/vcf/Filtered_VCFFILE.vcf.gz -e 'QUAL<=40' results/vcf/VCFFILE.vcf
```

### 11. Index filtered vcffiles

* Note: required installing tabix. Followed instructions here : https://github.com/trinityrnaseq/Griffithlab_rnaseq_tutorial_wiki/blob/master/AWS-Setup.md to install to SoilCocciSeq/results/vcf
* Set export path as:
export PATH=$PATH:/global/scratch/users/lcouper/SoilCocciSeqs/results/vcf/tabix-0.2.6

```
# in the SoilCocciSeq/results/vcf directory:
tabix *.vcf.gz # e.g. this should be the filtered vcf file
```





### 12. Create dictionary for reference genome

Software used: bio/picard/3.0.0-gcc-11.4.0
Command:

```
picard CreateSequenceDictionary -R RefGenome/CocciRef_GCA_000149335.2.fna
```

### 13. Sort vcf according to refeference dictionary 

Softwared used: bio/picard/3.0.0-gcc-11.4.0    
Script: sortvcf.sh   

```
picard SortVcf \
-I Filtered_VCFFILE.vcf.gz \
-O Filtered_Sorted_VCFFILE.vcf.gz \
-SD /global/scratch/users/lcouper/SoilCocciSeqs/RefGenome/CocciRef_GCA_000149335.2.dict
```


### 14. Filter low quality and rare SNPs using vcftools 

Software used: bio/vcftools/0.1.16-gcc-11.4.0, bio/bcftools/1.16-gcc-11.4.0    
Script name: filtersnps.sh    
Relevant snippet:   

```
vcftools --gzvcf Filtered_Sorted_VCFFILE.vcf.gz --maf 0.01 --minQ 30 --max-missing 0.90 --minDP 10 --recode --recode-INFO-all --out VCF_AllVariants.vcf
```
*kept 243,596 out of a possible 634596 Sites  
*or 156654 when run on just my 3 soil samples   
*or 9113 if you set -maf 0.05


#### 15. Remove multi-allelic sites

Software used: bio/vcftools/0.1.16-gcc-11.4.0, bio/bcftools/1.16-gcc-11.4.0   
Script name: biallelic.sh   
Relevant snippet:
```
bcftools view -m2 -M2 -v snps VCF_AllVariants.vcf > VCF_Biallelic.vcf 
```

To identify number of SNPs in vcf file:
grep -v "^#" VCF_Biallelic.vcf|wc -l
*Here, 232,417 SNPs remained*
*or 151843 when run on just my 3 sapmles 
*or 8930 for maf 0.05

### 16. Generate genotype matrix

Software used: bio/vcftools/0.1.16-gcc-11.4.0
Scriptname: genotypematx.sh
Relevant snippet: 

```
vcftools --012 --vcf Filtered_Sorted_VCFFILE_SJV_Genomes.vcf --out SJV_genotype_matr
```
*Note, this command outputs 3 files: ‘.012’ contains the genotypes of each individual on a separate line (with 0, 1, 2 denoting the number of non-reference alleles at the site), ‘.ind’ lists the individuals included in the main file, ‘.pos’ details the site location included in the main file.* 

#### LIC : 
for some reason, the genotype matrix output from above only included 1s and 0s... which is not correct.
So I then tried querying the VCFFILE and outputing this by doing the following: 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' VCF_Biallelic.vcf > output.txt
And playing with the resulting output.txt file in R...


### LIC:

look to this paper for next steps including PCA, fastADMIXTURE, using Nucmer to remove repetitive regions, etc:
https://journals.asm.org/doi/full/10.1128/mbio.01976-19 


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



