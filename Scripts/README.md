# 🧬 Processing *Coccidioides* Whole-Genome Sequences

![Pipeline](https://img.shields.io/badge/Pipeline-WGS-blue)
![Platform](https://img.shields.io/badge/Platform-HPC-green)
![Language](https://img.shields.io/badge/Scripts-Bash%20%7C%20R-orange)
![Tools](https://img.shields.io/badge/Tools-GATK%20%7C%20BWA%20%7C%20Samtools-purple)

This repository documents the scripts and steps used to process *Coccidioides* sequencing data from raw reads through downstream genomic analyses.

Note that the aligned genomes can be found under NCBI BioProject PRJNA1522484. 

---

## 🔬 Table of Contents

**[1. Raw data processing](#1-raw-data-processing)**

- [1.1 Obtain raw reads](#11-obtain-raw-reads)
- [1.2 Filter low-quality reads and trim bases](#12-filter-low-quality-reads-and-trim-bases)
- [1.3 Normalize read lengths to 75 bp](#13-normalize-read-lengths-to-75-bp)
- [1.4 Optional: Quality control with FastQC](#14-optional-quality-control-with-fastqc)

**[2. Alignment and BAM processing](#2-alignment-and-bam-processing)**

- [2.1 Align reads to reference genome](#21-align-reads-to-reference-genome)
- [2.2 Sort alignments and convert to BAM](#22-sort-alignments-and-convert-to-bam)
- [2.3 Optional: Extract mapping and coverage statistics](#23-optional-extract-mapping-and-coverage-statistics)
- [2.4 Add or replace read groups](#24-add-or-replace-read-groups)
- [2.5 Optional: Verify read groups and compute depth](#25-optional-verify-read-groups-and-compute-depth)
- [2.6 Mark and remove duplicates](#26-mark-and-remove-duplicates)
- [2.7 Index BAM files](#27-index-bam-files)

**[3. Variant calling](#3-variant-calling)**

- [3.1 Call variants using GATK HaplotypeCaller](#31-call-variants-using-gatk-haplotypecaller)
- [3.2 Combine GVCF files](#32-combine-gvcf-files)
- [3.3 Joint genotyping to produce metaVCF](#33-joint-genotyping-to-produce-metavcf)
- [3.4 Filter variants to produce project-specific VCF](#34-filter-variants-to-produce-project-specific-vcf)
- [3.5 Convert final vcf file to a pseudo-diploid genotype](#35-convert-final-vcf-file-to-a-pseudo-diploid-genotype)

**[4. Downstream genomic analyses](#4-downstream-genomic-analyses)**

- [4.1 Assess population structure: ADMIXTURE](#41-assess-population-structure-admixture)
- [4.2 Mating type locus assignment](#42-mating-type-locus-assignment)
- [4.3 Fst between environmental clusters](#43-fst-between-environmental-clusters)
- [4.4 Diversity metrics](#44-diversity-metrics)
- [4.5 Construct phylogenetic tree](#45-construct-phylogenetic-tree)
- [4.6 Assessing and correcting for linkage disequilibrium](#46-assessing-and-correcting-for-linkage-disequilibrium)
- [4.7 Twisst, window-based genomic relationships](#47-twisst-window-based-genomic-relationships)

---

**Software used**
- vcftools/0.1.16-gcc-11.4.0
- sratoolkit.3.0.7
- bwa-mem2/2.2.1
- samtools/1.17-gcc-11.4.0
- Trimmomatic V 0.39 (Bolger et al. 2014)
- fastp v 1.0.1 (manually installed from [here](https://github.com/OpenGene/fastp).
- gatk
- java
- python3
- picard/3.0.0-gcc-11.4.0
- fastqc/0.12.1-gcc-11.4.0

---

## 1. Raw data processing

#### 1.1 Obtain raw reads
**Notes:** Either download raw reads from sequencer or from NCBI SRA for previously published sequences.  
**Code:**   
```
# To download [this cocci genome](https://www.ncbi.nlm.nih.gov/sra/?term=SRR25635497):
~/Downloads/sratoolkit*/bin/prefetch SRR25635497

# To convert from .sra format to fastq:
~/Downloads/sratoolkit*/bin/fasterq-dump --split-files SRR25635497

# This will generate two files: filename_1.fastq and filename_2.fastq for paired-end reads
```

#### 1.2 Filter low-quality reads and trim bases
**Notes:** Illumina adapters [are available and downloaded from here (https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE.fa). Ensure this adapter sequence file is in the same folder as your fastq files     
**Code:**       
```
module load bio/trimmomatic/0.39-gcc-11.4.0
trimmomatic PE PS02PN14-1_S1_L007_R1_001.fastq.gz PS02PN14-1_S1_L007_R2_001.fastq.gz \
PS02PN14-1_S1_L007_R1_001.trim.fastq.gz PS02PN14-1_S1_L007_R1_001.untrim.fastq.gz \
PS02PN14-1_S1_L007_R2_001.trim.fastq.gz PS02PN14-1_S1_L007_R2_001.untrim.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:35 SLIDINGWINDOW:4:15
```

#### 1.3 Normalize read lengths to 75 bp 

**Notes:** This is to account for variation in sequenced read lengths across genomes (e.g. 150 bp, 300 bp, 75b p). We want to normalize to the lowest common denominator (here, 75 bp)  
**Code:**
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

#### 1.4 Optional: Quality control with FastQC

**Notes:** Outputs HTML report of per base and per sequence quality scores
**Code:**  
```
fastqc trimmed_fastqc/*.fastq.gz
```

---

## 2. Alignment and BAM processing

#### 2.1 Align reads to reference genome

**Notes:** Purpose is to determine where in the genome a given sequence/read is located. We used C. immitis RS [on NCBI here](https://www.ncbi.nlm.nih.gov/search/all/?term=GCA_000149335.2), after marking and removing duplicates.         
**Code:** 
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

#### 2.2 Sort alignments and convert to BAM

**Notes:** Compress sam to bam and sort bam file        
**Code:**
```
# Set $YOUR-PATH to point to where the gatk package lives
java -jar "$YOUR-PATH/gatk-package-4.5.0.0-local.jar" SortSam \
-I temp.sam \
-O temp.bam \
-SORT_ORDER coordinate
```

#### 2.3 Optional: Mapping and coverage statistics
**Notes:** Obtain mapping and coverage statistics for each sample    
**Code:**
```
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

#### 2.4 Add or replace read groups
**Notes:** The purpose of this step is to organize sequence data by library prep batch and sequencing runs parameters. Information on this step informed by these resources: [1](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)
[2](https://gatk.broadinstitute.org/hc/en-us/community/posts/4412745467931-HaplotypeCaller-does-not-work)    
**Code:** 
```
# Note this is for a single sample
picard AddOrReplaceReadGroups \
I=$YOUR-PATH-IN/PS02PN14-2_S2_L007.deduped.bam \
O=$YOUR-PATH-OUT/PS02PN14-2_S2_L007.rg.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=PS02PN14-2
```

#### 2.5 Optional: Compute depth
**Notes:** This step computes the depth at each position of the genome. This creates a txt file where the second and third columns are the position and coverage, respectively    
**Code:**  
```
# To compute depth at each position of the genome:   
samtools depth -a "$sample" > "${sample}.depth.txt"
done

# To calculate the mean depth from this file:
awk 'BEGIN { total = 0; count = 0 } { total += $3; count += 1; } END { avg = total / count; print avg} ' sample.depth.txt 
```

#### 2.6 Mark and remove duplicates
**Notes:** Duplicates reflect same sequence fragment being amplified and read multiple times. Keeping duplicates can lead to inflated estimates of coverage and can bias variant-calling steps    
**Code:**. 
```
picard MarkDuplicates \
-REMOVE_DUPLICATES TRUE \
-I ${sample}.sorted.bam \
-O "${sample}.deduped.bam" \
-M "${sample}.dup_metrics.txt"
done
```

#### 2.7 Index BAM files
**Notes:** Indexing creates faster referencing/extraction of specific regions based on genomic coordinates
**Code:**
```
samtools index ${sample}.deduped.bam
```

---

## 3. Variant calling

#### 3.1 Call variants using GATK HaplotypeCaller

**Notes:** Purpose is to identify specific positions in the genome of each sample that differ from the reference.    
**Code:**
```
java -jar "$YOUR-PATH/gatk-package-4.5.0.0-local.jar" HaplotypeCaller \
-R CocciRef_GCA_000149335.2.masked.fna \
-ploidy 1 \
-ERC GVCF \ 
-I 14B1.deduped.bam \
--output-mode EMIT_VARIANTS_ONLY \
-O 14B1.g.vcf.gz
```

#### 3.2 Combine GVCF files 
**Notes:** Purpose of this step is to create a dataset where all variant sites across all samples are considered.    
**Code:**
```
# first, combined all individual g.vcf.gz files created in the prior step into a single directory, then created a list of files in this directory:

ls *.vcf.gz > gvcfs.list

java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" CombineGVCFs \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant gvcfs.list \
-O combined.g.vcf.gz
```

#### 3.3 Joint genotyping to produce metaVCF 
**Notes:** Now, genotype samples using variant sites considered across all samples to create a global ('meta') vcf file, from which we can subset for specific analyses     
**Code:**
```
java -jar "$YOUR-PATH/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" GenotypeGVCFs \
-R CocciRef_GCA_000149335.2.masked.fna \
-ploidy 1 \
-V combined.g.vcf.gz \
-O metavcf.gz
```

#### 3.4 Filter variants to produce project-specific VCF
**Notes:** Subset the vcf to the samples included in a particular analyses. Then, flag and remove variants based on quality score, coverage, missingness etc. for just those samples    
**Code:**
```
# List all samples to include for a given analysis
Subset_analysis1.txt 

# Identify variants
java -jar "$YOUR-PATH/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" VariantFiltration \
-R CocciRef_GCA_000149335.2.masked.fna \
--variant metavcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || QUAL < 20" \
--filter-name "BasicAndBiasFilters" \
-O joint.filtered.vcf.gz

# Select (remove filtered) variants
java -jar "$YOUR-PATH//gatk-package-4.5.0.0-local.jar" SelectVariants \
-R CocciRef_GCA_000149335.2.masked.fna \
--variant joint.filtered.vcf.gz \
--restrict-alleles-to BIALLELIC \
--select-type-to-include SNP \
--exclude-filtered \
-O final.vcf.gz

# Keep only sites with >=90% samples genotyped 
vcftools --gzvcf final.vcf.gz \
  --max-missing 0.9 \
  --recode --recode-INFO-all \
  --out final_filtered_maxmissing
```


#### 3.5 Convert final vcf file to a pseudo-diploid genotype 
**Notes:** This step useful because haploid genotypes are not natively supported by vcftools and other packages    
**Code:** 
```
# First, create a 'ploidy' file to tell vcftools which part of the chromsome to consider haploid. Here, we are specificying all positions (by using large value of 999999999)
echo "* 0 999999999 . 2" > ploidy.txt

# Next, use the bcftools plug-in to correct ploidy across all sites (as specificed in the ploidy.txt file above)
bcftools +fixploidy final_filtered_maxmissing.vcf -- -p ploidy.txt > final.diploid.vcf
```

---

## 4. Downstream genomic analyses

### 4.1 ADMIXTURE
**Notes:** ADMIXTURE estimates ancestry to each of K clusters. We ran ADMIXTURE across K = 2-10 in our analysis. Note that the ADMIXTURE software requires that samples are in PLINK format. See [here](https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html) for more info   
**Code:**      
```
for K in 2 3 4 5 6 7 8 9 10; do
  for rep in $(seq 1 20); do
    seed=$((1000 + K * 100 + rep))
    run_prefix="K${K}_rep${rep}"
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

### 4.2 Mating type locus assignment 
**Notes:** Each isolate of *Coccidioides* has a mating type locus with one or two idiomorphs, MAT1-1 or MAT1-2, and sexual reproduction can only occur between distinct idiomorphs. Identifying the mating type locus for each individual and population can therefore provide clues about sexual reproduction and recombination. We first downloaded the MAT domain proteins from NCBI:
- [α-box domain (MAT1-1-1), C. immitis. EF472259.1](https://www.ncbi.nlm.nih.gov/search/all/?term=EF472259.1).
- [HMG domain (MAT1-2-1) from C. posadasii. EF472258.1](https://www.ncbi.nlm.nih.gov/search/all/?term=EF472258.1).
Then queried samples against these sequences. 
**Scripts:** See 'matingtype_updated.sbatch   


### 4.3 Fst between environmental clusters
**Notes:** Here we estimate genomic differentiation between different inferred clusters (using only the soil-derived isolates)   
**Code:**  
```
# First, note the samples to include in each 'population' in a txt file
echo -e "22AC2\n34B2\n58B1\n87A1" > Pop1_K23.txt
echo -e "13B1\nPS02PN14-1\n118a3\n157b2\nL100\n239a3b2" > Pop2_K23.txt

vcftools --vcf Subset_envr.final.diploid.vcf \
    --weir-fst-pop Pop1_K23.txt \
    --weir-fst-pop Pop2_K23.txt \
    --out fst_cc
```

### 4.4 Diversity metrics

**Notes:** Here we will calculate: the number of segregating sites (S), Nucleotide diversity (θπ), Watterson's theta (θw), and Tajima's D. note that because we are calculating metrics using only variant sites (ie from the VCF), we need to normalize based on the number of 'callable regions'. See **Script:** extract_callable_regions.py, which creates file: callable_regions.bed    
**Code:** 
```
# Number of segregating sites (S)
S_all=$(vcftools --vcf allsamples.final.recode.vcf --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
echo $S_all

# Nucleotide diversity (θπ). Average # of pairwise differences *per site* between all sequences in a population.
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

vcftools --vcf allsamples.final.diploid.vcf \
  --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_all_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_all_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_all =", s/c}'


# Watterson's theta (θw) (S, normalized by # of samples)
S=$(vcftools --vcf allsamples.final.recode.vcf --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(bcftools query -l allsamples.final.recode.vcf | grep -vxFf Clones.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"all theta_W: {(S / a_n) / callable}")
EOF

# Tajima's D
vcftools --vcf allsamples.final.diploid.vcf \
  --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_all
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_all =", sum/n}' tajimasD_all.Tajima.D
```

### 4.5 Construct phylogenetic tree 
**Notes:** We used the C. posadasii Silveira strain [SRR9644374](https://www.ncbi.nlm.nih.gov/biosample/?term=SRS007089) to root the tree   
**Code:**  
```
# First, convert vcf to phylipp:   
python3 ~/software/vcf2phylip.py \
  -i allsamples_withCpSilv.final.recode.vcf \
  -o CpSilv

# Next, test different models of molecular evolution   
iqtree3 -s allsamples_withCpSilv.final.recode.min4.phy \
        -m TESTONLY+ASC \ # test diff. nucleotide substituion models and pick the best one based on BIC. Includes ascertainment bias (For using VCF)
        -mset JC,HKY,K80,TN,GTR \
        -mrate E,G \
        -nt 8 \
        -o CpSilv \
        -pre modeltest_ASC_noR_rootCpSilv \
        -redo

# Finally, run tree using the best model as determined above
iqtree3 -s allsamples_withCpSilv.final.recode.min4.phy \
        -m GTR+F+ASC+G4 \ # best modeled determined in step 2 
        -bb 1000 \ #  1,000 bootstraps
        -alrt 1000 \ # 1,000 replicates of an approximate likelihood ratio test (to assess branch support)
        -nt 8 \ # run on 8 threads
        -o CpSilv \ # use C posadasii as outgroup
        -pre final_withCpSilv.bestmodel_1000 \
        -redo # overwrite old output
```

### 4.6 Assessing and correcting for linkage disequilibrium

**Notes:** We use the plink package the investigate LD and remove linked SNPs   
**Code:** 
```
# Convert filtered SNP VCF into PLINK binary format 
plink --vcf allsamples.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out allsamples_plink

# Create the LD-pruned dataset 
plink --bfile allsamples_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out allsamples_ld_r05

#  Make a pruned vcf file 
plink --bfile allsamples_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --recode vcf --out allsamples_ld_r05_pruned
```

### 4.7 Twisst, window-based genomic relationships
**Notes:** Twisst (topology weighting by iterative sampling of sub-trees) is a tool that assesses support for different genealogical relationships in moving windows along the genome. Note that twisst requires numpy and a tree-building software (we used raxml/8.2.12)   
**Code:**
```
# First, parse the vcf to .geno for twisst processing:
module load anaconda3/2024.10-1-11.4
export PYTHONPATH=$PYTHONPATH:$HOME/software/genomics_general
python $HOME/software/genomics_general/VCF_processing/parseVCF.py \
  -i cocci37.snps.vcf.gz \
  --ploidy 1 \
  --skipIndels \
  -o cocci37.geno.gz

# Then, create trees in sliding windows (this is the time-consuming step)   
python $HOME/software/genomics_general/phylo/raxml_sliding_windows.py \
  -g cocci37.geno.gz -p cocci37.w100 \
  --windType sites -w 100 -M 50 --model GTRCAT \  # 100 SNP sliding window
  --raxml raxmlHPC-AVX -T $SLURM_CPUS_PER_TASK \
  --log cocci37.w100.raxlog.txt

# Next, assess tree weights. Here, we will randomly sample 100 SNPs with replacement within each 100 SNP window. We'll do this for 100 bootstrap iterations. Then we'll repeat the twisst weighting using each iteration.  
See scripts linked below as these are longer code chunks
```
**Scripts:** boot_geno.py, twisst_summarise.py, twisst_boot.sbatch
