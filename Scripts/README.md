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
--filter-name "BasicAndBiasFilters" \
-O joint.filtered.vcf.gz
```

Step 2: "Select" (remove filtered) Variants
Software used: java, gatk 4.5.0.0     
Script: selectsnps.sbatch     
Code snippet:     

```
java -jar "/global/scratch/users/lcouper/SoilCocciSeqs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar" SelectVariants \
-R ../RefGenome/CocciRef_GCA_000149335.2.masked.fna \
--variant jointvcf_filtered.vcf.gz \
--restrict-alleles-to BIALLELIC \
--select-type-to-include SNP \
--exclude-filtered \
-O final.vcf.gz
```

### 16. Keep only sites with >=90% genotyped samples

Software used: vcftools/0.1.16-gcc-11.4.0
Script: NA (command line. Runs fast)
Code snippet:     
Note: Repeat for final_withCp.vcf

```
module load bio/vcftools/0.1.16-gcc-11.4.0
vcftools --gzvcf final.vcf.gz \
  --max-missing 0.9 \
  --recode --recode-INFO-all \
  --out final_filtered_maxmissing
```

Number of SNPs retained: 64,471. With C. posadasii : 64,767



### 17. Construct phylogenetic tree 

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

iqtree3 -s final_withCp.min4.phy \
        -m TEST \ # test different nucleotide substituion models and pick the best one based on BIC
        -bb 1000 \ #  1,000 bootstraps
        -alrt 1000 \ # 1,000 replicates of an approximate likelihood ratio test (to assess branch support)
        -nt AUTO \ # automatically detect and use number of optimal threads
        -o CpSilv \ # Specify C. posadasii Silveira strain as the outgroup
        -pre final_withCpSilv.min4_1000
```

One option for visualizing tree (but note we visualized in R using ggree):   
https://itol.embl.de/tree/136152214211185591747337347

### 18. Assess population structure 

Conducted using STRUCTURE v 2.3.4
Downloaded version for MacOS without front end [here](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html): 
On local computer, created directory 'structure_run' to store package and data files.  
To run structure analysis (must be in 'structure_run' directory.   
```
chmod +x structure  # May be necessary to run first if getting permission denied errors
./structure -m mainparams2 -K 2 -o output
```
Note the parameters listed in mainparams, and the format of the tester2.str are very specific. Copy/follow the versions attached here when running this analysis for real.    

Note that I did provide population identifiers, based on population structure observed in the PCA

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

# Additional downstream analyses 

### Scaffolding SNPs into genes 

Download .gtf file for cocci reference [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149335.2/).

### Identify size of each chromosome

```
cat CocciRef_GCA_000149335.2.fna | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```

![image](https://github.com/user-attachments/assets/3086e222-c492-4028-8700-0adc5b3c5ded)


## Examining mating type distribution

Step 1. Download gtf (and fna) files for each MAT idiomorphs from NCBI.   
For MAT1-1, I used the [C. immitis RS assemebly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000149335.2/).
For MAT1-2, I used the [C. immitis strain H538.4 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000149815.1/).

Step 2. Identify the genomic positions of the MAT locus for each representative
```
grep -i MAT-1 CimmitisRS_MAT1_1_Rep.gff > MAT_coords_mat1_1.gff
grep -i MAT1 CimmitisH538.4_MAT1_2_Rep.gtf > MAT_coords_mat1_2.gff
# Note the slightly different naming
```

Step 3. Index each assembly and extract the MAT loci using the (mRNA) coordinates identified above
```
samtools faidx Cimmitis_RS.fna
samtools faidx Cimmitis_H5384.fna
samtools faidx Cimmitis_RS.fna DS016985.1:384173-385395 > MAT1-1_RS.fna
samtools faidx Cimmitis_H5384.fna GG704913.1:1653763-1656710 > MAT1-2_H5384.fna
# Concatenate for mapping
cat MAT1-1_RS.fna MAT1-2_H5384.fna > MAT_combined.fna

# Note, for simplicity I fixed headers to: >MAT1_1_RS   >MAT1_2_H5384 (done manually in text editor)

# Index combined file
bwa-mem2 index MAT_combined.fna
```

Step 4. Identify most likely idiomorph for each cocci genome   
Here, we will align reads from each of the cocci genomes (ours and all others included in the analysis) to the MAT loci. Then we will identify to which idiomorph the coverage on the alignment is higher (indicating higher similarity).    
Software used: bwa-mem2/2.2.1, samtools/1.17-gcc-11.4.0       
Script: matingtype_loop.sbatch         
Note the output (the coverage to MAT1_RS or MAT_2_H5384) was directed to .tsv files (matingtype_coverage).    


Issue: results not matching up with [Engelthaler et al. 2016](https://journals.asm.org/doi/full/10.1128/mbio.00550-16#figS9). Potentially need to correct for mapping to other portions of these assemblies (eg rpb1)



### Alternative approach to mating type locus investigation 

Using just the protein sequences.
Downloaded the alpha-box protein sequence for MAT1-1 idiomorph (from C. immitis RS) [EF472259.1](https://www.ncbi.nlm.nih.gov/nuccore/EF472259.1).
and for MAT1-2 idiomorph (from C. posadasii) [EF472258.1](https://www.ncbi.nlm.nih.gov/nuccore/EF472258.1).
```
module load python/3.10.12-gcc-11.4.0 
module load spades/4.1.0 # Note that spades requires a more recent version of python
spades.py -1 SJV_9_1.fastq -2 SJV_9_2.fastq -o SJV_9_spades_output
```


## Fst differentiation between clinical and environmental isolates

First, created pop1 and pop2 txt files indicating assignment to environmental or clinical 'populations'. I first focused on only California samples to avoid spurious detection due to demographic processes. I then re-ran the analysis using only Washington samples to investigate how the Fst-outlier loci identified for California compared to those identified for Washington. 

California isolates:
```
echo -e "13B1\n14B1\n22AC2\n22BC1\34B2\n58B1\nPS02PN14-1\nPS02PN14-2\nPS02PN14-3" > CApop1.txt
echo -e "SD_1\nSJV_1\nSJV_10\nSJV_11\nSJV_2\nSJV_3\nSJV_4\nSJV_5\nSJV_6\nSJV_7\nSJV_8\nSJV_9\nUCLA293\nUCLA294\nUCLA295" > CApop2.txt
```
Then, I converted my final.vcf file to a pseudo-diploid genotype (as haploid genotypes are not natively supported by vcftools)
```
# First, create a 'ploidy' file to tell vcftools which part of the chromsome to consider haploid. Here, we are specificying all positions (by using large value of 999999999)
echo "* 0 999999999 . 2" > ploidy.txt

# Next, use the bcftools plug-in to correct ploidy across all sites (as specificed in the ploidy.txt file above)
module load bio/bcftools/1.16-gcc-11.4.0
bcftools +fixploidy final_filtered_maxmissing.recode.vcf -- -p ploidy.txt > final_diploid.vcf
```
Lastly, run vcftools to estimate Fst along the genome.   
Here, we estimated Fst for each SNP along the genome.   
Note: I previously tried estimating Fst in tiled windows, then calculating an empirical p-value using the permuted Fst distribution. But the issue was that I had to make the window size really large (>10 kbp) to get any significant hits after FDR correction. Given that, it seemed better to estimate per-SNP, and do a threshold-based comparision (i.e. candidate SNPs are those >99.9% CI from the permuted distribution.

```
cd /global/scratch/users/lcouper/SoilCocciSeqs/FinalOutputs
module load bio/vcftools/0.1.16-gcc-11.4.0
vcftools --vcf final_diploid.vcf \
    --weir-fst-pop CApop1.txt \
    --weir-fst-pop CApop2.txt \
    --fst-window-size 100000 \
    --out fst_100kbp_window_CA
```

Repeat for Washington isolates (may not keep as these environmental/clinical isolates are nearly clonal)
```
echo -e "A391\nA432\nA501\nA502\nWA_202\nWA_205\nWA_211\nWA_212\nWA_221" > WApop1.txt
echo -e "WA_1\nB11019\nB11034\nB12398\nB13956\nB15317\nB16692\nB17554" > WApop2.txt
vcftools --vcf final_diploid.vcf \
    --weir-fst-pop WApop1.txt \
    --weir-fst-pop WApop2.txt \
    --fst-window-size 100000 \
   --out fst_100kbp_window_WA
```

To assess statistical significance, randomly re-shuffle 'population' labels, and re-estimate Fst (repeat 500 times).    

```
cd /global/scratch/users/lcouper/SoilCocciSeqs/FinalOutputs

module load bio/vcftools/0.1.16-gcc-11.4.0

# number of permutations
nperm=500

# file with all sample IDs
samples=all_samples.txt

# number of samples in group 1 (e.g., environmental)
group1_n=10

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


## Tajima's D 

Tajima's D provides evidence of different types of selection. Its calculation is based on the site frequency spectrum.   

| Interpretation of Tajima's D | Description |
|-----------------------------|-------------|
| **Negative Tajima's D**     | Excess of rare alleles — may indicate purifying or positive selection. |
| **Positive Tajima's D**     | Excess of intermediate-frequency alleles — may indicate balancing selection. |
| **Tajima's D ≈ 0**          | Consistent with neutral evolution under constant population size. |   

   
Here, we want to calculate Tajima's D separately for the clinical and environmental samples. It is typically calculcated in  windows. I tried various window sizes but 100 kb seemed to be best 

Software used: vcftools/0.1.16-gcc-11.4.0
Code snippet (run at command line, very fast):
```
vcftools --vcf final_diploid.vcf \ # Note, requires this 'diploid' version as input
  --keep CApop2.txt \ # Names of the clinical CA samples stored in this text file
  --TajimaD 100000 \
  --out tajimasD_clinical

vcftools --vcf final_diploid.vcf \ # Note, requires this 'diploid' version as input
  --keep CApop1.txt \ # Names of the environmental CA samples stored in this text file
  --TajimaD 100000 \
  --out tajimasD_environmental
```

## Nucleotide diversity, θπ

θπ is the average number of pairwise differences *per site* between all sequences in a population.   
Here, we want to calculate θπ separately for clinical and environmental isolates from CA.    

| Interpretation of θπ | Description |
|-----------------------------|-------------|
| **High θπ in environmental isolates**     | Large, diverse populations. | 
| **Low θπ in clinical isolates**     | Selection / adpatation in the host; population bottlenecks. |    
| **Higher θπ in clinical vs environmental**     | Balancing selection in host |    

Software used:  vcftools/0.1.16-gcc-11.4.0     
Code snippet (run at command line, very fast):  
```
# For environmental isolates
vcftools \
  --vcf final_diploid.vcf \
  --keep CApop1.txt \
  --window-pi 100000 \
  --window-pi-step 100000 \
  --out pi_environmental

# For clinical isolates
vcftools \
  --vcf final_diploid.vcf \
  --keep CApop2.txt \
  --window-pi 100000 \
  --window-pi-step 100000 \
  --out pi_clinical
```

Note: We also calculated Fst, and θπ per SNP, rather than in tiled windows. Then we averaged per-SNP values by gene, using the genome annotation file for C. immitis RS (CimmitisRS.gtf). This part was done in R.


## pN/pS calculation 

- Requires multi-sample FASTA, reference genome (RefGenome/CocciRef_GCA_000149335.2.fna), and reference genome annotation file (RefGenome/genomic.gff)
- Note: this calculation is only done *within* a population (e.g. just the environmental isolates from a single population). In the code below, I've specified specific samples on which to run the calculation. Need to update this for the real version.

**Step 1. Extract coding sequence coordinates by gene**   
This uses the gene annotation file. Note: Make sure to MERGE overlapping intervals by gene (so just one set of coordinates per gene).
Script used: pnps/extract_and_merge_cds_coords.sh    
Code snippet:
```
#!/bin/bash

# Set input and output paths
GFF="RefGenome/genomic.gff"
OUT="cds_coords_merged.bed"

# Step 1: Extract CDS features
grep -P '\tCDS\t' "$GFF" > cds_features.gff

# Step 2: Convert to BED format (0-based start, 1-based end)
awk 'BEGIN{OFS="\t"} {
  split($9, a, /[;=]/);
  gene_id = "NA";
  for (i=1; i<=length(a); i++) {
    if (a[i] ~ /[Gg]ene[Ii][Dd]/) {
      gene_id = a[i+1];
      break
    } else if (a[i] ~ /^Parent$/) {
      gene_id = a[i+1];
    }
  }
  print $1, $4 - 1, $5, gene_id
}' cds_features.gff > cds_coords_raw.bed

# Step 3: Sort BED
sort -k1,1 -k2,2n cds_coords_raw.bed > cds_coords_sorted.bed

# Step 4: Merge overlapping CDS intervals per gene
bedtools merge -i cds_coords_sorted.bed -c 4 -o distinct > "$OUT"

# Clean up intermediate files (optional)
rm cds_features.gff cds_coords_raw.bed cds_coords_sorted.bed

echo "Done. Output written to $OUT"
```

**Step 2. Generate consensus genomes per sample**. 
For each sample, apply its variants (from a multisample VCF) to the reference genome to generate a personalized FASTA — i.e., the consensus genome.   
Script used: generate_consensus.sh   
Code snippet:
```
#!/bin/bash
#SBATCH --job-name=test_consensus
#SBATCH --account=fc_envids
#SBATCH --partition=savio3
#SBATCH --time=10:00:00
#SBATCH --output=test.out

## commands to run:

cd /global/scratch/users/lcouper/SoilCocciSeqs/

module load bio/bedtools2/2.31.0-gcc-11.4.0
module load bio/bcftools/1.16-gcc-11.4.0

VCF=FinalOutputs/final_filtered_maxmissing.vcf.gz
REF=RefGenome/CocciRef_GCA_000149335.2.masked.fna
SAMPLES=("13B1" "14B1")
OUTDIR=/global/scratch/users/lcouper/SoilCocciSeqs/pnps/consensus_genomes_test

mkdir -p "$OUTDIR"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Generating consensus genome for $SAMPLE..."
    bcftools view -c1 -s "$SAMPLE" -Oz -o "$OUTDIR/$SAMPLE.vcf.gz" "$VCF"
    bcftools index "$OUTDIR/$SAMPLE.vcf.gz"

    bcftools consensus -f "$REF" "$OUTDIR/$SAMPLE.vcf.gz" > "$OUTDIR/$SAMPLE.genome.fa"
done

echo "✅ Step 1 complete: consensus genomes in $OUTDIR"
```






