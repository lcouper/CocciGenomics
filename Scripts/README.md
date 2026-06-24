# 🧬 Processing *Coccidioides* Whole-Genome Sequences

![Pipeline](https://img.shields.io/badge/Pipeline-WGS-blue)
![Platform](https://img.shields.io/badge/Platform-HPC-green)
![Language](https://img.shields.io/badge/Scripts-Bash%20%7C%20R-orange)
![Tools](https://img.shields.io/badge/Tools-GATK%20%7C%20BWA%20%7C%20Samtools-purple)

This repository documents the scripts and steps used to process *Coccidioides* sequencing data from raw reads through downstream genomic analyses.

---

## 🔬 Table of Contents

## 1. Reference Genome Preparation (run once)

- [1.1 Mask repeats in reference genome](#11-mask-repeats-in-reference-genome)  
- [1.2 Index reference genome](#12-index-reference-genome)  

## 2. Raw Data Preparation

- [2.1 Obtain raw reads from Berkeley QB3](#21-obtain-raw-reads-from-berkeley-qb3)  
- [2.2 Download published sequences from NCBI SRA](#22-download-published-sequences-from-ncbi-sra)  
- [2.3 Filter low-quality reads and trim bases](#23-filter-low-quality-reads-and-trim-bases)  
- [2.4 Normalize read length to 75 bp](#24-normalize-read-length-to-75-bp)  
- [2.5 Optional: Quality control with FastQC](#25-optional-quality-control-with-fastqc)  

### 3. Alignment and BAM Processing

- [3.1 Align reads to reference genome](#31-align-reads-to-reference-genome)  
- [3.2 Sort alignments and convert to BAM](#32-sort-alignments-and-convert-to-bam)  
- [3.3 Optional: Extract mapping and coverage statistics](#33-optional-extract-mapping-and-coverage-statistics)  
- [3.4 Add or replace read groups](#34-add-or-replace-read-groups)  
- [3.5 Optional: Verify read groups and compute depth](#35-optional-verify-read-groups-and-compute-depth)  
- [3.6 Mark and remove duplicates](#36-mark-and-remove-duplicates)
- [3.7 Optional: Calculate genome coverage at >10× depth](#37-optional-calculate-genome-coverage-at-10-depth)  
- [3.8 Index BAM files](#38-index-bam-files)

### 4. Variant Calling

- [4.1 Call variants using GATK HaplotypeCaller](#41-call-variants-using-gatk-haplotypecaller)  
- [4.2 Combine GVCF files](#42-combine-gvcf-files)  
- [4.3 Joint genotyping to produce metaVCF](#43-joint-genotyping-to-produce-metavcf)  
- [4.4 Filter variants to produce project-specific VCF](#44-filter-variants-to-produce-project-specific-vcf)
- [4.5 Convert final vcf file to a pseudo-diploid genotype](#45-convert-final-vcf-file-to-a-pseudo-diploid-genotype)

### Downstream Genomic Analyses

- [5.1 FST differentiation between clinical and environmental isolates](#fst-differentiation-between-clinical-and-environmental-isolates)  
- [5.2 Population structure analysis](#assess-population-structure)  
- [5.3 Map SNPs to genes](#scaffolding-snps-into-genes)  
- [5.4 Determine chromosome sizes](#identify-size-of-each-chromosome)  
- [5.5 Mating type analysis](#examining-mating-type-distribution)  
- [5.6 Tajima’s D](#tajimas-d)  
- [5.7 Nucleotide diversity (θπ)](#nucleotide-diversity-θπ)  
- [5.8 McDonald–Kreitman test](#mk-test)  
- [5.9 Gene function and GO term analysis](#investigating-gene-function-and-go-terms)  
- [5.10 Extract amino acid sequences for differentiated genes](#get-amino-acid-sequence-for-significantly-differentiated-genes)  
- [5.11 Construct phylogenetic tree](#construct-phylogenetic-tree)
- [5.12 Linkage Disequilibrium](#linkage-disequilibrium)

---

**Software used**
- vcftools/0.1.16-gcc-11.4.0
- bio/bwa-mem2/2.2.1
- bio/samtools/1.17-gcc-11.4.0
- Trimmomatic V 0.39 (Bolger et al. 2014)
- fastp v 1.0.1 (manually installed from [here](https://github.com/OpenGene/fastp).
- gatk
- java
- python3
- bio/picard/3.0.0-gcc-11.4.0
- bio/fastqc/0.12.1-gcc-11.4.0


## Reference Genome Preparation
#### 1.1 Mask repeats in reference genome   
*Only need to do once*
Using reference genome for Coccidioides immitis RS (GCA_000149335.2)   
Downloaded [here for immitis](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/).     
and [here for posadasii](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000150055.1/).    
Using repeat library available [here](https://github.com/hyphaltip/cocci_repeats/).    
Purpose: Repetitive regions can lead to issues with reference alignment and variant calling    
Software used: repeatmasker/4.1.0   *Note: need to use this and not the newest version*   
Script: repeastmasker.sh *Note: this step done on SCG instead of Savio*    
Command:

```
RepeatMasker -pa 16 -lib immitis_repeats.fa --norna CocciRef_GCA_000149335.2.fna
```

#### 1.2 Index reference genome  
*Only need to do once*
Purpose: Enables quick access to specific locations of the genome (like the index of a book)   
   Script name: bwamem_index.sh      
Command:      

```
bwa-mem2 index CocciRef_GCA_000149335.2.fna.masked
```

Alternatively,  
```
samtools faidx CocciRef_GCA_000149335.2.masked.fna
```

## Raw data processing
#### 2.1 Obtain raw reads from Berkeley QB3

The fastq.gz files (1 forward, 1 reverse) are stored here on the Remais Group Shared Drive:
SPORE/WGS/Sequence data (All)/

and Berkeley's HPC BRC at:
`/global/scratch/users/lcouper/SoilCocciSeqs`

#### 2.2 Download published sequences from NCBI SRA

Prior Coccidioides sequences were downloaded from NCBI using the SRA toolkit.

Additional notes [here](https://docs.google.com/document/d/1gkM7m6TjQAOO1pwxe4X2DrIuMPuA3uGd6UalImpb-h4/edit).
Tracker for downloaded sequences and metadata [here](https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit).


#### 2.3 Filter low-quality reads and trim bases

Note that Illumina adapters [available and downloaded from here](https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE.fa). Ensure this adapter sequence file is in the same folder as your fastq files.   
Job Script: trim.sh and trim.sra.sh   
Relevant code snippet:       
```
module load bio/trimmomatic/0.39-gcc-11.4.0
trimmomatic PE PS02PN14-1_S1_L007_R1_001.fastq.gz PS02PN14-1_S1_L007_R2_001.fastq.gz \
PS02PN14-1_S1_L007_R1_001.trim.fastq.gz PS02PN14-1_S1_L007_R1_001.untrim.fastq.gz \
PS02PN14-1_S1_L007_R2_001.trim.fastq.gz PS02PN14-1_S1_L007_R2_001.untrim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```

#### 2.4 Normalize read lengths to 75 bp 

Note: this is because there is variation in sequenced read lengths across genomes (ours are all 150bp paired end, but prior genomes vary from 75 - 300 bp PE). We want to normalize to the lowest common denominator -- here 75 bp.

Scripts: run_fastp_len75.sbatch, run_fastp_len75_b.sbatch       
Relevant code snippet:
```
fastp \
  -i PS02PN14-1_S1_L007_R1_001.trim.fastq \
  -I PS02PN14-1_S1_L007_R2_001.trim.fastq \
  -o PS02PN14-1_S1_L007_R1_001.len75.trim.fastq \
  -O PS02PN14-1_S1_L007_R2_001.len75.trim.fastq \
  --max_len1 75 --max_len2 75 \
  --length_required 75 \
  --html PS02PN14_fastp_report.html --thread 4
```

#### 2.5 Optional: Quality control with FastQC

Script: fastqc.sh, fastqc.sra.sh        
Relevant code snippet:      
```
module load bio/fastqc/0.12.1-gcc-11.4.0
fastqc trimmed_fastqc/*.fastq.gz
```

## Alignment and BAM processing
#### 3.1 Align reads to reference genome

Purpose: To determine where in the genome a given sequence/read is located    
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

#### 3.2 Sort alignments and convert to BAM

Compress sam to bam and sort bam file     
Script name: SamToBam.sh, SamToBamSRA.sh

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SortSam \
-I temp.sam \
-O temp.bam \
-SORT_ORDER coordinate
```

#### 3.3 Optional: Extract mapping and coverage statistics

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

#### 3.4 Add or replace read groups

Followed guidance [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)
and issue was diagnosed [here](https://gatk.broadinstitute.org/hc/en-us/community/posts/4412745467931-HaplotypeCaller-does-not-work). See [this spreadsheet](https://docs.google.com/spreadsheets/d/1wrwSLeURp-E7LDD0SKT1wXEnrET5IziknmJWmXCB_7o/edit?gid=1963297784#gid=1963297784) for what read group parameters were added:

Purpose: Organize sequence data by library prep batch and sequencing runs parameters   
 Script name: addrg_loop.sbatch, addrgsrsa_loop.sbatch   (note to run in loop: upload tsv with read group info for all samples)   
Otherwise (single sample): addrg.sbatch, addrgsra.sbatch   
Code snippet for single sample:  

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

#### 3.5 Optional: Verify read groups and compute depth

To verify read groups added correctly:
```
samtools view -H results/bam/14B1.rg.bam 
```

To compute depth at each position of the genome:   
Software used: bio/samtools/1.17-gcc-11.4.0   
Script name: computedepth.sbatch, computedepth_sra.sbatch  
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

#### 3.6 Mark and remove duplicates

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

#### 3.7 Optional: Calculate genome coverage at >10× depth

Software used: bio/bedtools2/2.31.0-gcc-11.4.0, bio/samtools/1.17-gcc-11.4.0    
Script: depth10x.sbatch, depth10x_sra.sbatch    
Code snippet:    
```
# Note, we are only considering 'callable' bases in our count here
# i.e., excluding masked bases

for bam in "$BAM_DIR"/*.deduped.bam; do
  sample=$(basename "$bam" .deduped.bam)
  echo "Processing $sample..."

  # ensure BAM is indexed (skip if .bai exists)
  [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] || samtools index -@ 8 "$bam"

  # stream depths only within callable regions; count total callable bases and ≥10×
  percent=$(samtools depth -a -@ 8 -b "$CALLABLE_BED" "$bam" \
    | awk 'BEGIN{tot=0; ge10=0} {tot++; if($3>=10) ge10++} END{if(tot>0) printf("%.2f", 100*ge10/tot); else print "NA"}')

  echo -e "${sample}\t${percent}" >> "$OUTFILE"
done

```

*Note that any samples with <90% of their genome covered at >10x were then removed from downstream analysis. This is to avoid missing genotypes distorting population structure

#### 3.8 Index BAM files
Software used: bio/samtools/1.17-gcc-11.4.0   
Script: index_dedupedbams.sbatch, index_dedupedbams.sra.sbatch
Code snippet:
```
samtools index results/bam/${base}.deduped.bam
```

### Variant Calling

#### 4.1 Call variants using GATK HaplotypeCaller

Note on GATK installation: downloaded gatk from [here](https://github.com/broadinstitute/gatk/releases) and then uploaded the jar file to savio to working directory. Guidance on these steps found [here](https://www.biostars.org/p/405702/).   
Software used: java, gatk 4.5.0.0    
Script name: haplo.sh, haplosra.sh          
Code snippet:   

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" HaplotypeCaller \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
-ploidy 1 \
-ERC GVCF \ 
-I results/bam/58B1.deduped.bam \
--output-mode EMIT_VARIANTS_ONLY \
-O results/haplocalled/58B1.g.vcf.gz
```

#### 4.2 Combine GVCF files 

First, combined all the above files into a single directory 'AllGenomesHaploCalled'. Then, created a list of files in this directory using:
```
cd AllGenomesHaploCalled
ls *.vcf.gz > gvcfs.list
ls *.vcf.gz > gvcfs_withCp.list # repeat with CpSilv (outgroup for trees)

```
Purpose: Creates a dataset where all variant sites across all samples are considered. This enables variant callers to use information from one sample to infer the most likely genotype in another, improving sensitivity and accuracy in low coverage regions, and reducing false positives.
Here, all samples are included in the 'gvcfs.list'. We will filter the metaVCF later (step 15) for analyses using specific subsets of samples.  
Software used: java, gatk 4.5.0.0    
Script name: combinegvcfs.sbatch and combinegvcfs_WithCpSilv.sbatch (for phylo tree later)      
Note that this is done in batches of sample because otherwise the memory is exhausted. The original, not batched version of the script is: combinegvcfs_og.sbatch   
Code snippet:  

```
module load java

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" CombineGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant gvcfs.list \
-O combined.g.vcf.gz

```

#### 4.3 Joint genotyping to produce metaVCF

Software used: java, gatk 4.5.0.0   
Script name: genotypegvcfs.sh and genotypegvcfs_WithCpSilv.sh   

```
module load java
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" GenotypeGVCFs \
-R "/global/scratch/users/lcouper/SoilCocciSeqs/RefGenome/CocciRef_GCA_000149335.2.masked.fna" \
-ploidy 1 \
-V combined.g.vcf.gz \
-O metavcf.gz
```

#### 4.4 Filter variants to produce project-specific VCF

Here, we subset the vcf to the samples included in a particular analyses. Then, flag and remove variants based on quality score, coverage, missingness etc. for just those samples.

Subset_envrclin.txt contains the names for all environmental and (new) clinical samples (no repreps)
Subset_envrclin_Cp.txt contains the names for all environmental and (new) clinical samples (no repreps)
Subset_envr.txt contains the names for all environmental samples   (no repreps)
Subset_envr_withrepreps.txt contains the names for all environmental samples   (with repreps)
Subset_envrclinlegacy.txt contains the names for all environmental, (new) clinical samples, and legacy clinical samples (no repreps)
Subset_all_withCpSilv.txt contains the names for all samples (no preps) including Cp Silv


Scripts: 
- filtervcfs_clinenvr.sbatch (all of our environmental and clinical samples)
- filtervcfs_clinenvr_Cp.sbatch (our environmental and clinical samples, plus C. posadasii)
- filtervcfs_envr.sbatch (only our environmental samples)
- filtervcfs_envr_withrepreps.sbatch (only our environmental samples)
- filtervcfs_All.sbatch (all samples)
- filtervcfs_All_withCpSilv.sbatch (all samples, with Cp Silv as outgroup for phylo tree)

Software used: java, gatk 4.5.0.0, vcftools/0.1.16-gcc-11.4.0   
Code snippet:     

```
# Step 1: "Filter" (identify) Variants

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" VariantFiltration \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant jointvcf.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || QUAL < 20" \
--filter-name "BasicAndBiasFilters" \
-O joint.filtered.vcf.gz


# Step 2: "Select" (remove filtered) Variants

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SelectVariants \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant jointvcf_filtered.vcf.gz \
--restrict-alleles-to BIALLELIC \
--select-type-to-include SNP \
--exclude-filtered \
-O final.vcf.gz

# Step 3: Keep only sites with >=90% samples genotyped 

vcftools --gzvcf final.vcf.gz \
  --max-missing 0.9 \
  --recode --recode-INFO-all \
  --out final_filtered_maxmissing
```

Check how many SNPs retained:
```
# example:
bcftools view -H Subset_envrclin.final.recode.vcf | wc -l
```
metavcf: 263,266 
Subset_envr.final.recode.vcf: 62,847   
Subset_envr_withrepreps.final.recode.vcf: 63,377      
Subset_envrclin.final.recode.vcf: 56,791         
allsamples.final.recode.vcf: 56,201    
allsamples_withCpSilv.final.recode.vcf: 55,603   

#### 4.5 Convert final vcf file to a pseudo-diploid genotype 
Purpose: haploid genotypes are not natively supported by vcftools and other packages

```
# First, create a 'ploidy' file to tell vcftools which part of the chromsome to consider haploid. Here, we are specificying all positions (by using large value of 999999999)
echo "* 0 999999999 . 2" > ploidy.txt

# Next, use the bcftools plug-in to correct ploidy across all sites (as specificed in the ploidy.txt file above)
module load bio/bcftools/1.16-gcc-11.4.0
# run:
bcftools +fixploidy allsamples.final.recode.vcf -- -p ploidy.txt > allsamples.final.diploid.vcf
bcftools +fixploidy Subset_envr.final.recode.vcf -- -p ploidy.txt > Subset_envr.final.diploid.vcf
bcftools +fixploidy Subset_envrclin.final.recode.vcf -- -p ploidy.txt > Subset_envrclin.final.diploid.vcf
bcftools +fixploidy Subset_envrclin_Cp.final.recode.vcf -- -p ploidy.txt > Subset_envrclin_Cp.final.diploid.vcf
bcftools +fixploidy Subset_envr_withrepreps.final.recode.vcf -- -p ploidy.txt > Subset_envr_withrepreps.final.diploid.vcf
bcftools +fixploidy allsamples_withCpSilv.final.recode.vcf -- -p ploidy.txt > allsamples_withCpSilv.final.diploid.vcf
```

# Additional downstream analyses 

## Assess population structure: ADMIXTURE

Downlaoded ADMIXTURE [here](https://dalexander.github.io/admixture/download.html) and uploaded for use on savio  
Scripts: run_admixture_envrclin.sbatch, run_admixture_envrclin_Cp.sbatch (version with CpSilv)   
Code snippet:   
```
for K in 2 3 4 5 6 7 8 9 10; do
  for rep in $(seq 1 20); do
    seed=$((1000 + K * 100 + rep))
    run_prefix="K${K}_rep${rep}"
    echo "Running K=${K}, rep=${rep}, seed=${seed}"
    admixture --cv -s "$seed" -j8 "$admix_prefix.bed" "$K" | tee "${run_prefix}.log"
    mv "$(basename "$admix_prefix").${K}.Q" "${run_prefix}.Q"
    mv "$(basename "$admix_prefix").${K}.P" "${run_prefix}.P"
  done
done

echo -e "K\trep\tseed\tcv_error\tloglikelihood" > admixture_envrclin_replicate_summary.tsv

for log in K*_rep*.log; do
  K=$(echo "$log" | sed -E 's/K([0-9]+)_rep([0-9]+).log/\1/')
  rep=$(echo "$log" | sed -E 's/K([0-9]+)_rep([0-9]+).log/\2/')
  seed=$((1000 + K * 100 + rep))
  cv=$(grep "CV error" "$log" | awk '{print $4}')
  ll=$(grep "Loglikelihood" "$log" | tail -n 1 | awk '{print $2}')
  echo -e "${K}\t${rep}\t${seed}\t${cv}\t${ll}" >> admixture_envrclin_replicate_summary.tsv
done
```




## Assess population structure: STRUCUTRE 

Conducted using STRUCTURE v 2.3.4
Downloaded versions without front end [here](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html): 
Can run this on local computer or cluster (for multiple reps).  

On local computer: created directory 'structure_run' to store package and data files.  To run structure analysis, must be in 'structure_run' directory.   
```
chmod +x structure  # May be necessary to run first if getting permission denied errors
cd ~/Dropbox/CurrentProjects/SPORE/structure_run
./structure -m mainparams2 -K 2 -o output
```
Note the parameters listed in mainparams, and the format of the tester.str are very specific. Copy/follow the versions attached here when running this analysis for real.    

Note that I did provide population identifiers, based on population structure observed in the PCA, BUT I did not use this popinfo to inform the clustering in STRUCTURE (i.e. the clustering was unsupervised) 

Some relevant parameters from mainparams2:
```
Basic Program Parameters

#define MAXPOPS    2      // (int) number of populations assumed
#define BURNIN    10000   // (int) length of burnin period
#define NUMREPS   20000   // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   tester2.str   // (str) name of input data file
#define OUTFILE  outfile  //(str) name of output data file

Data file format

#define NUMINDS    76    // (int) number of diploid individuals in data file
#define NUMLOCI    196269    // (int) number of loci in data file
#define PLOIDY       1    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 1    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier
```

To test which level of K is most appropriate: 
```
# Reps 1–5
for rep in {1..5}; do ./structure -m mainparamsClinEnv -K 2 -D 2000${rep} -o output_clinenvr_K2_rep${rep}; done
for rep in {1..5}; do ./structure -m mainparamsClinEnv -K 3 -D 3000${rep} -o output_clinenvr_K3_rep${rep}; done
for rep in {1..5}; do ./structure -m mainparamsClinEnv -K 4 -D 4000${rep} -o output_clinenvr_K4_rep${rep}; done
for rep in {1..5}; do ./structure -m mainparamsClinEnv -K 5 -D 5000${rep} -o output_clinenvr_K5_rep${rep}; done
for rep in {1..5}; do ./structure -m mainparamsClinEnv -K 6 -D 6000${rep} -o output_clinenvr_K6_rep${rep}; done
```

To run on cluster (recommended).   
Script: run_structure_Kx_reps.sbatch (K = 1 through K = 8).   
With legacy genomes: run_structure_Kx_reps_legacy.sbatch (only conducted on K = 3 and K = 4).   


### Scaffolding SNPs into genes 

Download .gtf file for cocci reference [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/).

### Identify size of each chromosome

```
cat CocciRef_GCA_000149335.2.fna | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```
![image](https://github.com/user-attachments/assets/3086e222-c492-4028-8700-0adc5b3c5ded)


## Mating type locus assignment 

Each isolate of *Coccidioides* has a mating type locus with one or two idiomorphs, MAT1-1 or MAT1-2, and sexual reproduction can only occur between distinct idiomorphs. Identifying the mating type locus for each individual and population can therefore provide clues about sexual reproduction and recombination. 

Step 1. Download MAT domain proteins from NCBI (Note: downloaded on local computer, then uploaded to Savio)

[α-box domain (MAT1-1-1), C. immitis. EF472259.1](https://www.ncbi.nlm.nih.gov/search/all/?term=EF472259.1).

[HMG domain (MAT1-2-1) from C. posadasii. EF472258.1](https://www.ncbi.nlm.nih.gov/search/all/?term=EF472258.1).

Optionally, compare reuslts with [Engelthaler et al. 2016](https://journals.asm.org/doi/full/10.1128/mbio.00550-16#figS9) and [Teixeira et al. 2019](https://journals.asm.org/doi/full/10.1128/mbio.01976-19). 

Step 2. Query samples against these sequences.   
Script: matingtype_updated.sbatch   

## Fst differentiation between clinical and environmental isolates

First, created pop1a and pop1b txt files indicating assignment to environmental or clinical 'populations'. I focused on only samples with matching full ancestry (based on admixture results) to avoid spurious detection due to demographic processes. 

```
echo -e "22AC2\n22BC1\n34B2\n58B1\n87A1\n137a1_redo" > Pop1a.txt
echo -e "Kern6\nKern7\nKern13\nKern16\nKern21\nKern27" > Pop1b.txt
```

Then, run vcftools to estimate Fst along the genome.   
Here, we estimated Fst per site (can take averages by gene in R if desired)

```
vcftools --vcf Subset_envrclin.final.diploid.vcf \
    --weir-fst-pop Pop1a.txt \
    --weir-fst-pop Pop1b.txt \
    --out fst_per_site_Pop1
```

To assess statistical significance, randomly re-shuffle 'population' labels, and re-estimate Fst (repeat 500 times).     
Script: fst_perm.sbatch   
Code snippet:
```
# number of permutations
nperm=5000

# file with all sample IDs
samples=all_samples.txt

# number of samples in group 1 (e.g., environmental)
group1_n=15

mkdir -p perm_fst

for i in $(seq 1 $nperm); do
  echo "Permutation $i"

  # Shuffle and split samples
  shuf $samples > perm_fst/tmp_samples.txt
  head -n $group1_n perm_fst/tmp_samples.txt > perm_fst/group1.txt
  tail -n +$((group1_n + 1)) perm_fst/tmp_samples.txt > perm_fst/group2.txt

 # Run Fst
  vcftools --vcf final_diploid.vcf \
    --weir-fst-pop perm_fst/group1.txt \
    --weir-fst-pop perm_fst/group2.txt \
    --out perm_fst/fst_perm_$i \
    --stdout | grep -v "^#" | awk -v i=$i '{print $1, $2, $3, i}' >> perm_fst/fst_all_perms.txt
done
```

## Fst between environmental clusters

Similar to above, but using only environmental isolates and defining populations based on admixture output

```
echo -e "22AC2\n22BC1\n34B2\n58B1\n87A1\n137a1_redo" > Pop1.txt
echo -e "PS02PN14-1\nPS02PN14-2\nPS02PN14-3\n13B1\n14B1" > Pop2.txt
echo -e "118a3\n118b3\n157b2\n158b3\nL100\n239a3b2" > Pop3.txt
```

Then, run vcftools to estimate Fst along the genome.   
Here, we estimated Fst per site (can take averages by gene in R if desired)

```
vcftools --vcf Subset_envr.final.diploid.vcf \
    --weir-fst-pop Pop1.txt \
    --weir-fst-pop Pop2.txt \
    --out fstPop12

vcftools --vcf Subset_envr.final.diploid.vcf \
    --weir-fst-pop Pop1.txt \
    --weir-fst-pop Pop3.txt \
    --out fstPop13

vcftools --vcf Subset_envr.final.diploid.vcf \
    --weir-fst-pop Pop2.txt \
    --weir-fst-pop Pop3.txt \
    --out fstPop23
```
Results: 
Fst 1 & 2: Weir and Cockerham mean Fst estimate: 0.27331; weighted Fst estimate: 0.41696   
Fst 1 & 3: Weir and Cockerham mean Fst estimate: 0.15339; weighted Fst estimate: 0.22844
Fst 2 & 3: Weir and Cockerham mean Fst estimate: 0.22700; weighted Fst estimate: 0.31262

### Diversity metrics

These calculations will use the filtered VCF file that contains all samples (rather than subset-specific filtered vcf files since that will confounded diversity calculations due to QC steps).   
First, create txt files indicating sample names for each subset:   

```
echo -e "22AC2\n22BC1\n34B2\n58B1\n87A1\n137a1_redo\nPS02PN14-1\nPS02PN14-2\nPS02PN14-3\n13B1\n14B1\n118a3\n118b3\n157b2\n158b3\nL100\n239a3b2" > Envr.txt      
bcftools query -l Subset_envrclin.final.diploid.vcf | grep '^Kern' > Clin.txt   
cat Envr.txt Clin.txt > EnvrClin.txt
``` 

**S (number of segregating sites)**

```
S_envr=$(vcftools --vcf allsamples.final.recode.vcf --keep Envr.txt --mac 1 --recode --stdout | grep -vc "^#")
S_clin=$(vcftools --vcf allsamples.final.recode.vcf --keep Clin.txt --mac 1 --recode --stdout | grep -vc "^#")
S_envrclin=$(vcftools --vcf allsamples.final.recode.vcf --keep EnvrClin.txt --mac 1 --recode --stdout | grep -vc "^#")
S_all=$(vcftools --vcf allsamples.final.recode.vcf --mac 1 --recode --stdout | grep -vc "^#")

echo "S_envr = $S_envr"
echo "S_clin = $S_clin"
echo "S_envrclin = $S_envrclin"
echo "S_all = $S_all"
```
S_envr: 43,545\
S_clin: 51,191\
S_envrclin: 53,882\
S_all: 56,201   

**Watterson's theta (S, normalized by # of samples)**

```
# Environmental
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Envr.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(wc -l < Envr.txt)
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

python3 - <<EOF
S = $S
n = $n
callable = $callable
a_n = sum(1/i for i in range(1, n))
theta_w = (S / a_n) / callable
print(f"environmental theta_W: {theta_w}")
EOF

# Clinical
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Clin.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(wc -l < Clin.txt)
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

python3 - <<EOF
S = $S
n = $n
callable = $callable
a_n = sum(1/i for i in range(1, n))
theta_w = (S / a_n) / callable
print(f"clinical theta_W: {theta_w}")
EOF

# Environmental and clinical
S=$(vcftools --vcf allsamples.final.recode.vcf --keep EnvrClin.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(wc -l < EnvrClin.txt)
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

python3 - <<EOF
S = $S
n = $n
callable = $callable
a_n = sum(1/i for i in range(1, n))
theta_w = (S / a_n) / callable
print(f"envr and clin theta_W: {theta_w}")
EOF

# All (including legacies)
S=$(vcftools --vcf allsamples.final.recode.vcf --mac 1 --recode --stdout | grep -vc "^#")
n=$(bcftools query -l allsamples.final.recode.vcf | wc -l)
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

python3 - <<EOF
S = $S
n = $n
callable = $callable
a_n = sum(1/i for i in range(1, n))
theta_w = (S / a_n) / callable
print(f"all theta_W: {theta_w}")
EOF
```
environmental theta_W: 0.0005455180932584511
clinical theta_W: 0.0005681607521474858
envr and clin theta_W: 0.0005274288114176259
all theta_W: 0.00047252173477471156

**Nucleotide diversity, θπ**  
θπ is the average number of pairwise differences *per site* between all sequences in a population. **Key note: because we are calculating pi using only variant sites (ie from the VCF), we need to normalize based on the number of 'callable regions'.   
We did this using:extract_callable_regions.py (python script in RefGenonme directory) to create a file: callable_regions.bed. This calculation requires using diploid version of vcf.  
```
# environmental:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Envr.txt \
  --site-pi \
  --max-missing 0.9 \
  --out pi_environmental_sitewise

callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_environmental_sitewise.sites.pi)

awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_environmental =", s/c}'

# clinical
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Clin.txt \
  --site-pi \
  --max-missing 0.9 \
  --out pi_clinical_sitewise

callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_clinical_sitewise.sites.pi)

awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_clinical =", s/c}'

# envr and clin:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep EnvrClin.txt \
  --site-pi \
  --max-missing 0.9 \
  --out pi_envrclin_sitewise

callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_envrclin_sitewise.sites.pi)

awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_envrclin =", s/c}'

# all:
vcftools --vcf allsamples.final.diploid.vcf \
  --site-pi \
  --max-missing 0.9 \
  --out pi_all_sitewise

callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_all_sitewise.sites.pi)

awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_all =", s/c}'
```

pi_environmental = 0.000612615
pi_clinical = 0.000653024
pi_envrclin = 0.000666852
pi_all = 0.000674365


**Tajima's D**

Typically calculcated in  windows. I tried various window sizes but 100 kb seemed to be best. This calculation requires using diploid version of vcf. 
```
# Environmental
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Envr.txt \
  --TajimaD 100000 \
  --out tajimasD_environmental

awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_environmental =", sum/n}' tajimasD_environmental.Tajima.D


# Clinical
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Clin.txt \
  --TajimaD 100000 \
  --out tajimasD_clinical

awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_clinical =", sum/n}' tajimasD_clinical.Tajima.D


# Environmental + clinical
vcftools --vcf allsamples.final.diploid.vcf \
  --keep EnvrClin.txt \
  --TajimaD 100000 \
  --out tajimasD_envrclin

awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_envrclin =", sum/n}' tajimasD_envrclin.Tajima.D


# All samples
vcftools --vcf allsamples.final.diploid.vcf \
  --TajimaD 100000 \
  --out tajimasD_all

awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_all =", sum/n}' tajimasD_all.Tajima.D
```
mean_TajimasD_environmental = 1.24342  
mean_TajimasD_clinical = 1.17419  
mean_TajimasD_envrclin = 1.3271  
mean_TajimasD_all = 1.65487   


### MK Test

- Requires multi-sample FASTA, reference genome (RefGenome/CocciRef_GCA_000149335.2.fna), and reference genome annotation file (RefGenome/genomic.gff)
- Note: this calculation was done for only California genomes to minimize confounding by demography

**Step 1. (only need to run once) Extract coding sequence coordinates by gene from the reference**     
This uses the gene annotation file.     
Script used: pnps/extract_cds_bed12.sh. Original version: pnps/extract_and_merge_cds_coords.sh    
To run:
```
bash pnps/extract_cds_bed12.sh RefGenome/genomic.gff pnps/cds_coords_merged.bed12
```
**Step 2. Generate consensus genomes per sample** 
For each sample, apply its variants (from a multisample VCF) to the reference genome to generate a personalized FASTA — i.e., the consensus genome.   
Script used: generate_consensus_genomes.sh

**Step 3. Extract CDS sequences from each sample's consensus genome**

Software used: bedtools 2.31.0, bcftools 1.16     
Script: extract_sample_cds_from_consensus.sh. (Prior version: generate_per_sample_gene_vcfs_og.sh).   
Note that ~2% ambiguity is expected due to earlier repeat masking, and will be excluded in downstream analyses.   

**Step 4. Translate nucleotide sequences to proteins**  

Software used: biopython, python     
Script: translate_all_cds.py  (old version: fasta_to_protein.py)
*Run as: python translate_all_cds.py*

**Step 5. Filter problematic genes**   

Here, we are removing CDS with any of the following: internal stop codons (likely due to sequencing errors), >10% missing (coded as Xs), or lack of representation in >=75% of each group (clinical or environmental). We are removing those here as they will cause issues in downstream steps and/or introduce biases in our analyses. The output is one FASTA per gene containing only the passing samples.    
Note this first requires making a file that indicates to which group each sample belongs, called 'sample_to_group.tsv'  
Python script used: filter_genes.py or (without the 75% rule): filter_genes_no75rule.py  
*Run as: python filter_genes.py*

**Step 6. Align protein sequence per gene across samples**

Software used: muscle v3.8. Note the latest versions (v5) was giving issues, hence going with an older release. Program (muscle3.8.31_i86linux32.tar) was manually downloaded [here](https://drive5.com/muscle/downloads_v3.htm).   
Script used: run_muscle_alignments.py      
*Run as: run_muscle_alignments.py*    Note this may take ~30 minutes to run

**Step 7. Create per-gene CDS FASTA files across samples**    
I.e. we need the nucleotide sequence of each gene from each sample. This is required input for PAL2NAL.    
Script used: generate_cds_by_gene_strict.py
*Run as: generate_cds_by_gene_strict.py*   
Note there is a 'non-strict' version: generate_cds_by_gene.py that may hold on to duplicates.   

**Step 8. Create codon-aware nucleotide alignments**     
Software used: pal2nal      
Note I manually downloaded PAL2NAL from [here](https://www.bork.embl.de/pal2nal/#Download). Then uploaded to BRC, unpacked, and added to my path.       
Script used: run_pal2nal.sh   

**Step 9. Build codon level group consensuses (i.e. clinical vs environmental)**
Software used: python/3.10.12-gcc-11.4.0
Note this uses the 'sample_to_group.tsv' previously created, which lists samples and group assignments
Script uesd: build_group_consensus.py

**Step 10. Count within-group polymorphisms (pS/pS)**
Software used: python/3.10.12-gcc-11.4.0, egglib v3
Script used:  compute_pnps_by_group.py 
Note: may take ~10 minutes to run

**Step 11. Count between group "divergences" (dN/dS)**
Software used: python/3.10.12-gcc-11.4.0
Script used: compute_dnds_between_groups.py   

**Step 12. Conduct MK test from pN/pS and dN/dS counts**
Software used: python/3.10.12-gcc-11.4.0
Script used: mk_test_from_counts.py




### Investigating gene function and GO terms

Search Gene ID here on NCBI gene search, i.e. here: https://www.ncbi.nlm.nih.gov/gene/?term=Coccidioides+immitis+CIMG_02011 

### Get amino acid sequence for significantly differentiated genes

Purpose: Genes identified as significant on the basis of Fst could be the result of neutral or selective processes. Identifying whether there is corresponding amino acid changes at these genes can help resolve this.   
Script: aminoacid_pull.sbatch   # But note this only includes clinical samples -- need to add in environmental into the BAM call   


### Construct phylogenetic tree 

** Note: In order to root the phylogenetic tree, we used the C. posadasii Silveira strain [SRR9644374](https://www.ncbi.nlm.nih.gov/biosample/?term=SRS007089) **
These fastqs were then taken through the same steps as all other samples above (e.g. starting from step 1)   
The resulting vcf files were only used for rooting the tree.   
The vcf file *without* C. posadasii were used in all other analyses 


Step 1. Convert vcf to phylipp:   
Run at command line, very fast.  
Note that the 'vcf2phylip.py' script was downloaded from [here](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) and must be in working directory for command to work
```
python3 vcf2phylip.py -i allsamples_withCpSilv.final.recode.vcf -o CocciSamplesTree
# or for LD-pruned version:
python3 vcf2phylip.py -i allsampleswithCpSilv_ld_r05_pruned.vcf -o CocciSamplesTreeLD

```

Step 2. Test different models of molecular evolution   
Software used: iqtree/3.0.0    
Script: phylo_tree_testmodels.sh, phylo_tree_testmodels_LD.sh    
Code snippet:
```
iqtree3 -s allsamples_withCpSilv.final.recode.min4.phy \
        -m TESTONLY+ASC \ # test diff. nucleotide substituion models and pick the best one based on BIC. Includes ascertainment bias (For using VCF)
        -mset JC,HKY,K80,TN,GTR \
        -mrate E,G \
        -nt 8 \
        -o CpSilv \
        -pre modeltest_ASC_noR_rootCpSilv \
        -redo
```

Step 3. Run tree using the best model as determined in step 2    
Software used: iqtree/3.0.0    
Script: phylo_tree_bestmodel.sh    
Code snippet:
```
iqtree3 -s allsamples_withCpSilv.final.recode.min4.phy \
        -m GTR+F+ASC+G4 \ # best modeled determined in step 2 
        -bb 1000 \ #  1,000 bootstraps
        -alrt 1000 \ # 1,000 replicates of an approximate likelihood ratio test (to assess branch support)
        -nt 8 \ # run on 8 threads
        -o CpSilv \ # use C posadasii as outgroup
        -pre final_withCpSilv.bestmodel_1000 \
        -redo # overwrite old output
```

One option for visualizing tree (but note we visualized in R using ggree):   
https://itol.embl.de/tree/136152214211185591747337347



## Assessing and correcting for Linkage Disequilibrium ##

Downloaded plink version 1.9 [here](https://www.cog-genomics.org/plink/). (specifically the 64-bit Linux, stable beta version)
Uploaded folder to cluster and made it executable. 

Step 1. Convert filtered SNP VCF into PLINK binary format (can run at command line. very fast)    
```
cd /global/scratch/users/lcouper/SoilCocciSeqs/FinalOutputs
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --vcf allsamples.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out allsamples_plink

# with CpSilv (outgroup for tree)
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --vcf allsamples_withCpSilv.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out allsampleswithCpSilv_plink

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --vcf Subset_envrclin.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out Subset_envrclin_plink

# our environmental and clinical isolates plus Cp
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --vcf Subset_envrclin_Cp.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out Subset_envrclin_Cp_plink
```

Step 2. Create the LD-pruned dataset (again, can run at command line. very fast).    
Here we are using a window size of 50 SNPs, sliding by 5 SNPs each time. Within each window, plink identifies SNP pairs with r2 > 0.5 and removes variants until no remaining pair exceeds this threshold.
```
cd /global/scratch/users/lcouper/SoilCocciSeqs/FinalOutputs
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsamples_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out allsamples_ld_r05

# with CpSilv
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsampleswithCpSilv_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out allsampleswithCpSilv_ld_r05

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out Subset_envrclin_ld_r05

# our environmental and clinical isolates plus Cp
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_Cp_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out Subset_envrclin_Cp_ld_r05
```
This removed ~43,121 out of 56,201 variants, leaving 13,080 SNPs 

Step 3. Make pruned plink files (for any downstream analyses that use plink)
```
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsamples_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --make-bed --out allsamples_ld_r05_pruned

# with Cp Silv
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsampleswithCpSilv_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --make-bed --out allsampleswithCpSilv_ld_r05_pruned

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --make-bed --out Subset_envrclin_ld_r05_pruned

# our environmental and clinical isolates plus Cp
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_Cp_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --make-bed --out Subset_envrclin_Cp_ld_r05_pruned
```

Step 4. Make a pruned vcf file 
```
# with Cp Silv
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsampleswithCpSilv_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --recode vcf --out allsampleswithCpSilv_ld_r05_pruned

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --recode vcf --out Subset_envrclin_ld_r05_pruned
```
Note the pruned vcf is called 'allsamples_ld_r05_pruned.vcf' or 'Subset_envrclin_ld_r05_pruned.vcf'





