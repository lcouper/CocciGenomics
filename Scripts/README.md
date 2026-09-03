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
  - [4.4.1 Number of segregating sites (S)](#441-number-of-segregating-sites-s)
  - [4.4.2 Watterson's theta (θw)](#442-wattersons-theta-θw)
  - [4.4.3 Nucleotide diversity (θπ)](#443-nucleotide-diversity-θπ)
  - [4.4.4 Tajima's D](#444-tajimas-d)
- [4.5 Construct phylogenetic tree](#45-construct-phylogenetic-tree)
- [4.6 Assessing and correcting for linkage disequilibrium](#46-assessing-and-correcting-for-linkage-disequilibrium)
- [4.7 Twisst, window-based genomic relationships](#47-twisst-window-based-genomic-relationships)
- [4.8 Identifying deletion](#48-identifying-deletion)

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

These calculations will use the filtered VCF file that contains all samples (rather than subset-specific filtered vcf files since that will confound diversity calculations due to QC steps).   
First, create txt files indicating sample names for each subset:   

```
echo -e "22AC2\n22BC1\n34B2\n58B1\n87A1\n137a1_redo\nPS02PN14-1\nPS02PN14-2\nPS02PN14-3\n13B1\n14B1\n118a3\n118b3\n157b2\n158b3\nL100\n239a3b2" > Envr.txt      
bcftools query -l Subset_envrclin.final.diploid.vcf | grep '^Kern' > Clin.txt   
cat Envr.txt Clin.txt > EnvrClin.txt
```

For these diversity calculations, we exclude one isoalte from each pair that appear nearly clonal (i.e. <200 SNP differences)

| Near-clonal lineage | Keep | Drop (raw VCF name) | Pairwise SNP distance |
|---|---|---|---|
| Carrizo_1 / Carrizo_2 | 13B1 | **14B1** | 14 |
| PS02PN14 trio | PS02PN14-1 | **PS02PN14-2, PS02PN14-3** | 88 / 92 / 121 (3a–3b / 3c–3b / 3a–3c) |
| Bakersfield_1 / _2 | 22AC2 | **22BC1** | 25 |
| Coalinga_1 / _2 | 157b2 | **158b3** | 83 |
| McKittrick_1 / _2 | 118a3 | **118b3** | 120 |
| VFI19 / VFI20 | Kern19 | **Kern20** | 29 |
| VFI25 / VFI5 | Kern25 | **Kern5** | 79 |

```
# list out the near-clones to remove
cat > Clones.txt <<'EOF'
14B1
PS02PN14-2
PS02PN14-3
22BC1
158b3
118b3
Kern5
Kern20
EOF
```

#### 4.4.1 Number of segregating sites (S)

```
S_envr=$(vcftools --vcf allsamples.final.recode.vcf --keep Envr.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
S_clin=$(vcftools --vcf allsamples.final.recode.vcf --keep Clin.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
S_envrclin=$(vcftools --vcf allsamples.final.recode.vcf --keep EnvrClin.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
S_all=$(vcftools --vcf allsamples.final.recode.vcf --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")

S_envr_pop1=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop1.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
S_envr_pop2=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop2.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
S_envr_pop3=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop3.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")

# read out with: echo $S_envr
```


#### 4.4.2 Watterson's theta (θw)
*S, normalized by # of samples*

```
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

# Environmental   
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Envr.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt Envr.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"environmental theta_W: {(S / a_n) / callable}")
EOF

# Clinical    
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Clin.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt Clin.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"clinical theta_W: {(S / a_n) / callable}")
EOF

# Environmental and clinical   
S=$(vcftools --vcf allsamples.final.recode.vcf --keep EnvrClin.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt EnvrClin.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"envr and clin theta_W: {(S / a_n) / callable}")
EOF

# All (including legacies)    
S=$(vcftools --vcf allsamples.final.recode.vcf --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(bcftools query -l allsamples.final.recode.vcf | grep -vxFf Clones.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"all theta_W: {(S / a_n) / callable}")
EOF


# Soil population 1     
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop1.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt Pop1.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"soil pop1 theta_W: {(S / a_n) / callable}")
EOF

# Soil population 2     
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop2.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt Pop2.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"soil pop2 theta_W: {(S / a_n) / callable}")
EOF

# Soil population 3     
S=$(vcftools --vcf allsamples.final.recode.vcf --keep Pop3.txt --remove Clones.txt --mac 1 --recode --stdout | grep -vc "^#")
n=$(grep -vxFf Clones.txt Pop3.txt | wc -l)
python3 - <<EOF
S = $S; n = $n; callable = $callable
a_n = sum(1/i for i in range(1, n))
print(f"soil pop3 theta_W: {(S / a_n) / callable}")
EOF
```


#### 4.4.3 Nucleotide diversity (θπ)
θπ is the average number of pairwise differences *per site* between all sequences in a population. **Key note: because we are calculating pi using only variant sites (ie from the VCF), we need to normalize based on the number of 'callable regions'.   
We did this using:extract_callable_regions.py (python script in RefGenonme directory) to create a file: callable_regions.bed. This calculation requires using diploid version of vcf.  
```
callable=$(awk '{sum += $3 - $2} END {print sum}' ../RefGenome/callable_regions.bed)

# environmental:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Envr.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_environmental_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_environmental_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_environmental =", s/c}'

# clinical:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Clin.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_clinical_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_clinical_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_clinical =", s/c}'

# envr and clin:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep EnvrClin.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_envrclin_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_envrclin_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_envrclin =", s/c}'

# all (including legacies):
vcftools --vcf allsamples.final.diploid.vcf \
  --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_all_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_all_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_all =", s/c}'

# soil population 1:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop1.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_pop1_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_pop1_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_pop1 =", s/c}'

# soil population 2:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop2.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_pop2_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_pop2_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_pop2 =", s/c}'

# soil population 3:
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop3.txt --remove Clones.txt \
  --site-pi --max-missing 0.9 \
  --out pi_pop3_sitewise
sum_pi=$(awk 'NR > 1 {sum += $3} END {print sum}' pi_pop3_sitewise.sites.pi)
awk -v s="$sum_pi" -v c="$callable" 'BEGIN {print "pi_pop3 =", s/c}'
```

pi_environmental = 0.00063078   
pi_clinical = 0.000659787  
pi_envrclin = 0.000672379   
pi_all = 0.000674751   

pi_pop1 = 0.000461357   
pi_pop2 = 0.00032205    
pi_pop3 = 0.00052863    



#### 4.4.4 Tajima's D

Typically calculcated in  windows. I tried various window sizes but 100 kb seemed to be best. This calculation requires using diploid version of vcf. 
```
# Environmental
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Envr.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_environmental
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_environmental =", sum/n}' tajimasD_environmental.Tajima.D

# Clinical
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Clin.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_clinical
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_clinical =", sum/n}' tajimasD_clinical.Tajima.D

# Environmental + clinical
vcftools --vcf allsamples.final.diploid.vcf \
  --keep EnvrClin.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_envrclin
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_envrclin =", sum/n}' tajimasD_envrclin.Tajima.D

# All samples (including legacies)
vcftools --vcf allsamples.final.diploid.vcf \
  --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_all
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_all =", sum/n}' tajimasD_all.Tajima.D

# Soil population 1
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop1.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_pop1
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_pop1 =", sum/n}' tajimasD_pop1.Tajima.D

# Soil population 2
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop2.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_pop2
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_pop2 =", sum/n}' tajimasD_pop2.Tajima.D

# Soil population 3
vcftools --vcf allsamples.final.diploid.vcf \
  --keep Pop3.txt --remove Clones.txt \
  --TajimaD 100000 \
  --out tajimasD_pop3
awk 'NR > 1 && $4 != "nan" {sum += $4; n++} END {print "mean_TajimasD_pop3 =", sum/n}' tajimasD_pop3.Tajima.D
```

mean_TajimasD_environmental = 0.825641   
mean_TajimasD_clinical = 1.13959   
mean_TajimasD_envrclin = 1.19482    
mean_TajimasD_all = 1.60022    

mean_TajimasD_pop1 = 1.14279   
mean_TajimasD_pop2 = 2.25243    
mean_TajimasD_pop3 = 1.15089    


### 4.5 Construct phylogenetic tree 

** Note: In order to root the phylogenetic tree, we used the C. posadasii Silveira strain [SRR9644374](https://www.ncbi.nlm.nih.gov/biosample/?term=SRS007089) **
These fastqs were then taken through the same steps as all other samples above (e.g. starting from step 1)   
The resulting vcf files were only used for rooting the tree.   
The vcf file *without* C. posadasii were used in all other analyses 


Step 1. Convert vcf to phylipp:   
Run at command line, very fast.  
Note that the 'vcf2phylip.py' script was downloaded from [here](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) and added to:
/global/home/users/lcouper/Software
```
python3 ~/software/vcf2phylip.py \
  -i allsamples_withCpSilv.final.recode.vcf \
  -o CpSilv

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



### 4.6 Assessing and correcting for linkage disequilibrium

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

# just environmental isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --vcf Subset_envr.final.diploid.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out Subset_envr_plink
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

# just environmental
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envr_plink --allow-extra-chr --indep-pairwise 50 5 0.5 --out Subset_envr_ld_r05
```
This removed 44,564 out of 56,282 variants, leaving 11,718 SNPs 

Step 3. Make pruned plink files (for any downstream analyses that use plink)
```
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsamples_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --make-bed --out allsamples_ld_r05_pruned

# with Cp Silv
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsampleswithCpSilv_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --make-bed --out allsampleswithCpSilv_ld_r05_pruned

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --make-bed --out Subset_envrclin_ld_r05_pruned

# our environmental and clinical isolates plus Cp
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_Cp_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --make-bed --out Subset_envrclin_Cp_ld_r05_pruned

# just environmental
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envr_plink --allow-extra-chr --extract Subset_envr_ld_r05.prune.in --make-bed --out Subset_envr_ld_r05_pruned
```

Step 4. Make a pruned vcf file 
```
# with Cp Silv
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsampleswithCpSilv_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --recode vcf --out allsampleswithCpSilv_ld_r05_pruned

# with just our environmental and clinical isolates
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envrclin_plink --allow-extra-chr --extract Subset_envrclin_ld_r05.prune.in --recode vcf --out Subset_envrclin_ld_r05_pruned

# just envirionmental
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile Subset_envr_plink --allow-extra-chr --extract Subset_envr_ld_r05.prune.in --recode vcf --out Subset_envr_ld_r05_pruned

# all samples
/global/scratch/users/lcouper/SoilCocciSeqs/plink/plink --bfile allsamples_plink --allow-extra-chr --extract allsamples_ld_r05.prune.in --recode vcf --out allsamples_ld_r05_pruned
```
Note the pruned vcf is called 'allsamples_ld_r05_pruned.vcf' or 'Subset_envrclin_ld_r05_pruned.vcf'


### 4.7 Twisst, window-based genomic relationships

We are primarily interseted in our novel clinical samples, but we will use a few legacy clinical genomes as a control here. Therefore, we will use the vcf file with all samples (no re-preps) and CpSilv as an outgroup. These sample names live in: FinalOutputs/samples43.txt. We then subset the meta-sample vcf to these 43 using (which keeps only sites that are polymorphic within these 46 samples): 

```
bcftools view -S samples37.txt -m2 -M2 -v snps allsamples_withCpSilv.final.recode.vcf -Ou \
  | bcftools +fill-tags -Ou -- -t AC,AN \
  | bcftools view -e 'AC==0 || AC==AN' -Oz -o cocci37.snps.vcf.gz
``` 
Note that twisst requires numpy and a tree-building software (we will use raxml/8.2.12 here)

```
# from within the 'FinalOutputs/twisst' directory, parse the vcf to .geno for twisst processing (fast):
module load anaconda3/2024.10-1-11.4
export PYTHONPATH=$PYTHONPATH:$HOME/software/genomics_general
python $HOME/software/genomics_general/VCF_processing/parseVCF.py \
  -i ../cocci37.snps.vcf.gz \
  --ploidy 1 \
  --skipIndels \
  -o cocci37.geno.gz
```

Then, create trees in sliding windows:      
Script: twisst_trees.sbatch  
Code snippet:   
```
python $HOME/software/genomics_general/phylo/raxml_sliding_windows.py \
  -g cocci37.geno.gz -p cocci37.w100 \
  --windType sites -w 100 -M 50 --model GTRCAT \  # 100 SNP sliding window
  --raxml raxmlHPC-AVX -T $SLURM_CPUS_PER_TASK \
  --log cocci37.w100.raxlog.txt
```

Next, create the 'groupings' file. Note, these are the 'pure' soil isolates, that we'll use the 'reference groups' here, in addition to the Cp 'reference group'. 
```
mkdir -p groups out

# Reference groups are clone-corrected: one isolate per burrow, so leave-one-out is
# automatically leave-one-lineage-out. Patient clones are retained (they are focals, never
# references). Dropped: 22BC1, 14B1, 158b3, 118b3, PS02PN14-2, PS02PN14-3.

mkdir -p groups out

# Primary spec, K = 2: ADMIXTURE's Bakersfield vs everything else, clone-corrected.
cat > groups/refs_K2.tsv <<'EOF'
22AC2 Bakersfield
34B2 Bakersfield
58B1 Bakersfield
87A1 Bakersfield
13B1 NonBakersfield
PS02PN14-1 NonBakersfield
118a3 NonBakersfield
157b2 NonBakersfield
L100 NonBakersfield
239a3b2 NonBakersfield
CpSilv Outgroup
EOF

# Alternative spec: three groups defined by SITE rather than by ADMIXTURE. Justified by
# pairwise divergence -- Coalinga-Tracy (14,776) and Carrizo-McKittrick (15,926) are the two
# most similar site pairs, Bakersfield is most divergent from all (20,743-22,383).
cat > groups/refs_site3.tsv <<'EOF'
22AC2 Bakersfield
34B2 Bakersfield
58B1 Bakersfield
87A1 Bakersfield
13B1 CarrizoMcKittrick
PS02PN14-1 CarrizoMcKittrick
118a3 CarrizoMcKittrick
157b2 CoalingaTracy
L100 CoalingaTracy
239a3b2 CoalingaTracy
CpSilv Outgroup
EOF
```

Now we will run the tree generation step. Here, we will randomly sample 100 SNPs with replacement within each 100 SNP window. We'll do this for 100 bootstrap iterations. Then we'll repeat the twisst weighting using each iteration.  

First, created: boot_geno.py and twisst_summarise.py (see attached scripts. this lives in the twisst directory), and Scripts/twisst_boot.sbatch
Then, in the FinalOutputs/twisst directory: 
```
mkdir -p logs

for S in K2a K2b K3a K3b; do
  sbatch --export=ALL,SPEC=$S --array=1-100%25 ../../Scripts/twisst_boot.sbatch
done
```
After this finishes, download FinalOutputs/Twisst/boot_[K2a or K2b or K3a or K3b]/summaries to local machine (i.e. SPORE/Twisst/BootstrapSummaries.  



### 4.8 Identifying deletion

Three clinical isolates identified as unusual based on phylogenetic tree (they form a monophyletic clade that is basal to all other *C. immitis* and most similar to *C. posadasii*). To identify the underlying cause, we used the genotype matrix (GenotypeMatrix_ClinEnvr.csv) to look for sites where all three isolates were simultaneously `NA` and found 894 such sites. These NAs can mean either "no reads" (deletion) or that reads are too diverged to align. These were distinguished by using the bams for these three samples on Savio (stored in results/bam/Normalized75).    
Here, coverage *shape* is the informative statistic: a deletion gives a hard-edged block of zero depth, whereas divergence gives patchy coverage with conserved islands and graded boundaries. We included a non-deletion carrier in the same commands for comparison:

```
# mean depth in 500 bp bins across the candidate tract
samtools depth -a -r "$DEL1" "$bam" | awk '{b=int(($2-3154575)/500); s[b]+=$3; n[b]++}
  END {for(i=0;i<=102;i++) if(n[i]) printf "%6d %6.1f %s\n", 3154575+i*500, s[i]/n[i], (s[i]/n[i]<2?"ZERO":"")}'
```

We then confirmed the deletion and pinned its breakpoints based on discordant read pairs. Namely, if the region is deleted, read pairs spanning it should align with an insert size approximately equal to the deletion length. The filter must bracket the expected size — a first attempt using `$9>45000` for a ~39 kb deletion removed exactly the informative pairs.

```
samtools view "$bam" GG704911.1:3162000-3204000 | awk '$9>30000 && $9<50000'
```

This returned pairs with left mate ~3,163,0xx–3,163,2xx, right mate 3,202,77x–3,202,85x, insert ~39.6–39.8 kb, MAPQ 60 and many `NM:i:0`, along with 18–22% soft-clipping in the trio versus 3–6% in controls. Breakpoints: last covered base 3,163,754, first covered base after 3,202,775. The exact edge can also be read directly:
```
samtools depth -a -r GG704911.1:3162500-3164000 "$bam" | awk '$3>0{last=$2} END{print last}'   # expect 3163754
```

We found that the three isolates share two large deletions (~47 kb total) that are private to the three of them, and that this sequence is present in *C. posadasii*, so it is a derived loss rather than sequence they never had.
Coordinates reused throughout:
```
DEL1="GG704911.1:3163755-3202774"        # 39,020 bp, removes 6 protein-coding genes
D1_INTERIOR="GG704911.1:3170000-3200000" # for depth screens; avoids the breakpoints
DEL2="GG704913.1:4053986-4065345"        # ~8 kb deleted, with a retained ~1.75 kb island at 4,060,236-4,061,986
FLANK="GG704911.1:3300000-3351404"       # normalizer: same size as DEL1, similar masking
```

We then screened all other isolates for this deletion. Note that depth must be normalized to each sample's **own flank**, never compared raw across samples: the reference used for alignment is repeat-masked and these tracts are 36–42% `N` versus a 19–24% scaffold average, so even non-carriers sit at a ratio of ~0.4–0.6 rather than 1.0. Printing the flank value alongside the ratio also catches failed libraries, which otherwise produce false deletion calls. 

```
meandepth () { samtools depth -a -r "$2" "$1" | awk '{t+=$3;n++} END{if(n) printf "%.1f", t/n; else printf "0"}'; }

for bam in *.len75.aligned.sorted.deduped.bam; do
  s=${bam%%_*}; f=$(meandepth "$bam" "$FLANK"); a=$(meandepth "$bam" "$D1_INTERIOR")
  awk -v s="$s" -v f="$f" -v a="$a" 'BEGIN{r=(f>0?a/f:0); printf "%-14s flank=%7.1f del1=%7.1f (%.2f) %s\n", s,f,a,r,(r<0.10?"<== DEL1":"")}'
done | sort
```

Result: ratio = 0.00 for all three trio isolates and 0.27–0.63 for every other isolate. Both deletions are unique to the clade.   
Lastly, we sought to determine the putative function of genes inside the deletion. We downloaded the C. immitis RS annotation. We then overlapped the deletion interval with the GFF, finding 6 protein-coding genes, all annotated `hypothetical protein`: `CIMG_02019, CIMG_02016, CIMG_02011, CIMG_02010, CIMG_12729, CIMG_12730`. The neighboring CMGC/SRPK protein kinase `CIMG_02020` ends 84 bp before the breakpoint and is **intact**. Deletion 2 is gene-free in the current annotation. We then pulled the protein sequences (note that `protein.faa` headers key on `protein_id`, not `locus_tag`, and that `grep -A` on a FASTA ignores record boundaries)
```
zcat GCA_000149335.2_ASM14933v2_protein.faa.gz \
| awk '/^>/{keep = /EAS36666|EAS36665|EAS36662|EAS36657|EAS36656|KJF60083|KJF60084/} keep' > deleted_proteins.faa
```
Lastly, we identified whether this deletion was present in C. posadasii. To do so, we chopped the C. posadasii reference regions into 500 bp chunks and mapped them to the *C. posadasii* assembly. The flank serves as a positive control that the aligner crosses the *immitis*–*posadasii* boundary. Important: use the **unmasked** reference here. The masked version (used for read alignment) is ~42% `N` in these tracts, which silently destroys the test.

```
module load bio/samtools bio/bwa-mem2/2.2.1
REF=../../../RefGenome/CocciRef_GCA_000149335.2.fna
CP=../../../RefGenomeCp/GCF_018416015.2_ASM1841601v2_genomic.fna

chunk () { : > "$4"; for s in $(seq "$2" 500 "$3"); do e=$((s+499)); [ "$e" -gt "$3" ] && e="$3";
           samtools faidx "$REF" "$1:$s-$e" >> "$4"; done; }

chunk GG704911.1 3300000 3351404 flank_chunks.fa
bwa-mem2 mem -t 8 "$CP" flank_chunks.fa 2> log | samtools view -c -F 4 -F 256 -F 2048 -
```

| region | chunks mapping to *C. posadasii* |
|---|---|
| flank `GG704911.1:3,300,000-3,351,404` (control) | 103/103 = 100% |
| DEL1 `GG704911.1:3,163,755-3,202,774` | 50/79 = 63.3% |
| DEL2 `GG704913.1:4,053,986-4,065,345` | 23/23 = 100% |

We found that *C. posadasii* has the sequence, so it was present in the common ancestor and the trio's absence is a derived loss. 

To pull out coverage from this region to generate figure in R: 
```
module load bio/samtools
# Do this in the results/bam/Normalized75 folder with all the deduplicated .bam files and their index files
> deletion_depth.tsv
for bam in *.len75.aligned.sorted.deduped.bam; do
  case "$bam" in *reprep*) continue ;; esac
  s=$(echo "$bam" | sed -E 's/\.len75\..*//; s/_S[0-9]+_L[0-9]+$//')
  for reg in GG704911.1:3158755-3207774 GG704913.1:4048986-4070345; do
    samtools depth -a -r "$reg" "$bam" \
      | awk -v s="$s" 'BEGIN{OFS="\t"} {print s,$1,$2,$3}'
  done
done >> deletion_depth.tsv
gzip -f deletion_depth.tsv
```
