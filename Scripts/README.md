## Steps and scripts used to process _Coccidioides_ genomes from environmental samples or Sequence Read Archive ### 

### 1. Obtained raw reads from Berkeley QB3 or SRA.  

The fastq.gz files (1 forward, 1 reverse) are stored here on the Remais Group Shared Drive: SPORE/WGS/Sequence data (All)/ 
and Berkeley's HPC BRC at: /global/scratch/users/lcouper/SoilCocciSeqs

Prior cocci sequences were downloaded from NCBI using the SRA toolkit.   
Additional notes [here](https://docs.google.com/document/d/1gkM7m6TjQAOO1pwxe4X2DrIuMPuA3uGd6UalImpb-h4/edit?tab=t.0).   
Tracker for sequences downloaded and metadata [here](https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit?gid=1963297784#gid=1963297784). 

### 2. Filter poor quality reads and trim poor quality bases 

Note that Illumina adapters [available and downlaoded from here](https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE.fa). Ensure this adapter sequence file is in the same folder as your fastq files.   
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

### 4. Mask repeats in reference genome   

Using reference genome for Coccidioides immitis RS (GCA_000149335.2)   
Downloaded here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/  
Using repeat library available here:   
Purpose: Repetitive regions can lead to issues with reference alignment and variant calling    
Software used: repeatmasker/4.1.0   *Note: need to use this and not the newest version*   
Script: repeastmasker.sh *Note: this step done on SCG instead of Savio*    
Command:

```
RepeatMasker -pa 16 -lib immitis_repeats.fa --norna CocciRef_GCA_000149335.2.fna
```

### 4. Index reference genome  

Purpose: Enables quick access to specific locations of the genome (like the index of a book)   
Software used: bio/bwa-mem2/2.2.1    
Script name: bwamem_index.sh      
Command:      

```
bwa-mem2 index CocciRef_GCA_000149335.2.fna.masked
```

Alternatively, 
Software used: bio/samtools/1.17-gcc-11.4.0  
```
samtools faidx CocciRef_GCA_000149335.2.masked.fna
```

### 5. Align sequences to reference genome    

Purpose: To determine where in the genome a given sequence/read is located    
Software used: bio/bwa-mem2/2.2.1   
Script name: alignreads.sh, alignreads.sra.sh    
Relevant code snippet:   

```
#First unzip trimmed fastq files if not done already, then align to ref genome
gunzip trimmed_fastq/*.gz

for infile in trimmed_fastq/*_R1_001.trim.fastq
do
base=$(basename ${infile} _R1_001.trim.fastq)
bwa-mem2 mem -t 12 RefGenome/CocciRef_GCA_000149335.2.masked.fna \
trimmed_fastq/${base}_R1_001.trim.fastq trimmed_fastq/${base}_R2_001.trim.fastq > results/sam/"${base}.aligned.sam"
done
```

### 6. Sort and convert to bam

Compress sam to bam and sort bam file 
Software used: gatk
Script name: SamToBam.sh, SamToBamSRA.sh

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SortSam \
-I temp.sam \
-O temp.bam \
-SORT_ORDER coordinate
```

### 7. Extract mapping and coverage statistics 

Software used: samtools, java
Script name: MappingStats.sh, MappingStatsSRA.sh
Note that the cap for coverage is 250 so there may be a peak in the histograms at this value

```
# Loop through each sorted BAM file
for bam_file in "$bam_dir"/*.sorted.bam; do
    sample_id=$(basename "$bam_file" .sorted.bam)

# Mapping stats
samtools flagstat "$bam_file" > "$stats_dir/${sample_id}_mapping.stats.log"

# Coverage stats
picard CollectWgsMetrics \
  I="$bam_file" \
  O="$stats_dir/${sample_id}_coverage_stats.txt" \
  R="$ref"
done
```

### 8. Add or replace read groups 

Followed guidance [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)
and issue was diagnosed [here](https://gatk.broadinstitute.org/hc/en-us/community/posts/4412745467931-HaplotypeCaller-does-not-work). See [this spreadsheet](https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit?gid=1963297784#gid=1963297784) for what read group parameters were added:

Purpose: Organize sequence data by library prep batch and sequencing runs parameters   
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
RGSM=PS02PN14-2
```
### 9. Optional: Verify read groups and compute death 

To verify read groups added correctly:
```
samtools view -H results/bam/14B1.rg.bam 
```

To compute depth at each position of the genome:   
Software used: bio/samtools/1.17-gcc-11.4.0   
Script name: computedepth.sbatch
```
module load bio/samtools/1.17-gcc-11.4.0

for infile in *.aligned.sorted.bam
do
  echo "working with file $infile"
  base=$(basename "$infile" .aligned.sorted.bam)
  samtools depth -a "$infile" > "${base}.depth.txt"
done
```

Note: This creates a txt file where the second and third columns are the position and coverage, respectively.
To calculate the mean depth from this file:

```
awk 'BEGIN { total = 0; count = 0 } { total += $3; count += 1; } END { avg = total / count; print avg} ' results/bam/58B1.depth.txt
```


### 10. Mark and remove duplicates 

Purpose: Duplicates reflect same sequence fragment being amplified and read multiple times. Keeping duplicates can lead to inflated estimates of coverage and can bias variant-calling steps    
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


### 11. Index bam files with read group added

Software used: bio/samtools/1.17-gcc-11.4.0   
Script: index_dedupedbams.sbatch, index_dedupedbams.sra.sbatch
Code snippet:
```
samtools index results/bam/${base}.deduped.bam
done
```

### 12. Call variants using GATK HaplotypeCaller 

Note on GATK installation: downloaded gatk from [here](https://github.com/broadinstitute/gatk/releases) and then uploaded the jar file to savio to working directory. Guidance on these steps found [here](https://www.biostars.org/p/405702/).   
Software used: java, gatk 4.5.0.0    
Script name: haplo.sh, haplosra.sh          
Code snippet:   

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" HaplotypeCaller \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
-ploidy 1 \
-ERC GVCF \  # Note that this option specifies we only want SNPs retained
-I results/bam/58B1.deduped.bam \
--output-mode EMIT_VARIANTS_ONLY \
-O results/haplocalled/58B1.g.vcf.gz
```

### 13. Combine GVCF files 

First, combined all the above files into a single directory 'AllGenomesHaploCalled'. Then, created a list of files in this directory using:
```
cd AllGenomesHaploCalled
ls *.vcf.gz > gvcfs.list
```
Purpose: Creates a dataset where all variant sites across all samples are considered. This enables variant callers to use information from one sample to infer the most likely genotype in another, improving sensitivity and accuracy in low coverage regions, and reducing false positives   
Software used: java, gatk 4.5.0.0    
Script name: combinegvcfs.sh

```
module load java

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" CombineGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant gvcfs.list \
-O combined.g.vcf.gz

```

### 14. Joint-genotyping on combined GVCF files 

Software used: java, gatk 4.5.0.0   
Script name: genotypegvcfs.sh    

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" GenotypeGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
-ploidy 1 \
-V combined.g.vcf.gz \
-O jointvcf.vcf.gz
```

### 15. Flag and remove variants based on quality score, coverage, etc.

Step 1: "Filter" (identify) Variants
Software used: java, gatk 4.5.0.0     
Script: filtervcfs.sbatch     
Code snippet:     

```
module load java

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" VariantFiltration \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant jointvcf.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || QUAL < 20" \
--filter-name "BasicAndBiasFilters" \s
-O joint.filtered.vcf.gz
```

Step 2: "Select" (remove filtered) Variants
Software used: java, gatk 4.5.0.0     
Script: selectsnps.sbatch     
Code snippet:     

```
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SelectVariants \
-R ../RefGenome/CocciRef_GCA_000149335.2.fna \
--variant joint.vcf.filtered.vcf.gz \
--restrict-alleles-to BIALLELIC \ # most downstream techniques (GWAS, pop gen analyses) require biallelic SNP data
--select-type-to-include SNP \
-O final.SNPs.vcf.gz
```

Optional: unzip final.SNPs.vcf file and identify number of variant sites:

```
gunzip final.SNPs.vcf.gz
grep -v "^#" AllGenomesHaploCalled/final.SNPs.vcf | wc -l   # 196,269. With C. posadasii : 484,074
```


### 16. Construct phylogenetic tree 

** Note: In order to root the phylogenetic tree, we used the C. posadasii Silveira strain [SRR9644374](https://www.ncbi.nlm.nih.gov/biosample/?term=SRS007089) **
These fastqs were then taken through the same steps as all other samples above (e.g. starting from step 1)   
The resulting vcf files were only used for rooting the tree.   
The vcf file *without* C. posadasii were used in all other analyses 


Step 1. Convert vcf to phylipp:   
Run at command line, very fast.  
Note that the 'vcf2phylip.py' script was downloaded from [here](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) and must be in working directory for command to work
```
python3 vcf2phylip.py -i final_withCp.vcf -o CocciSamples
```

Step 2. Build phylogenetic tree    
Software used: iqtree/3.0.0    
Script: phylo_tree.sh    
Code snippet:
```
module load iqtree/3.0.0
# iqtree3 -s final.SNPs.min4.phy -m GTR+G -nt AUTO -o CpSilv (# non-bootstrapped, fast version)

iqtree3 -s final.SNPs.min4.phy \
        -m TEST \ # test different nucleotide substituion models and pick the best one based on BIC
        -bb 1000 \ #  1,000 bootstraps
        -alrt 1000 \ # 1,000 replicates of an approximate likelihood ratio test (to assess branch support)
        -nt AUTO \ # automatically detect and use number of optimal threads
        -o CpSilv \ # Specify C. posadasii Silveira strain as the outgroup
        -pre final.SNPs.withCpSilv.min4_1000
```

One option for visualizing tree (but note we visualized in R using ggree):   
https://itol.embl.de/tree/136152214211185591747337347

### 17. Assess population structure 

Conducted using STRUCTURE v 2.3.4
Downloaded version for MacOS without front end [here](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html): 
On local computer, created directory 'structure_run' to store package and data files.  
To run structure analysis (must be in 'structure_run' directory.   
```
./structure -m mainparams -K 2 -o output
```
Note the parameters listed in mainparams, and the format of the tester.str are very specific. Copy/follow the versions attached here when running this analysis for real.  



## Other steps: 

### Investigate SNP position in genes 

Download .gtf file for cocci reference [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/).

### Identify size of each chromosome

```
cat CocciRef_GCA_000149335.2.fna | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```

![image](https://github.com/user-attachments/assets/3086e222-c492-4028-8700-0adc5b3c5ded)







