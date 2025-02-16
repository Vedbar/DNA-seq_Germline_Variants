# Germline Variants Pipeline

This repository provides a step-by-step pipeline for processing germline variant analysis using genomic sequencing data. The pipeline includes raw data preprocessing, alignment, variant calling, and annotation using widely used bioinformatics tools.

---

## Table of Contents

1. [Server Information](#1-server-information)  
2. [Genomic Sequencing Data](#2-genomic-sequencing-data)
3. [QC raw reads](#3-qc-raw-reads)
4. [Align reads](#4-align-reads)
5. [Sort alignment](#5-sort-alignment)
6. [Mark duplicates](#6-mark-duplicates)
7. [Index BAM file](#7-index-bam-file)
8. [Variant Calling](#8-variant-calling)
9. [Combine Variants](#9-combine-variants)
10. [Genotype Variants](#10-genotype-variants)
11. [Variant Filtering](#11-variant-filtering)
12. [Genotype Concordance](#12-genotype-concordance)
13. [Variant Annotation](#13-variant-annotation)
14. [Extract Variant Fields](#14-extract-variant-fields)
15. [Class Exercise](#15-class-exercise)

---


## 1. Server Information

### Prerequisites
- **Tools:** Install  [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html), [FileZilla](https://filezilla-project.org/download.php), [IGV](https://igv.org/doc/desktop/)  
- **Server Details:**  
  - **IP Address:** `168.105.161.70`  
  - **Port:** `22`  
  - **Access:** Requires JABSOM or UH network (use VPN for remote access).
  - Note: With PuTTY and FileZilla you can connect to server.
 

### Security Practices
- Avoid multiple failed login attempts to prevent account locking.  
- Use strong passwords or SSH keys.  
- Log out after completing tasks.

### Installing Miniconda
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Installing software
```bash
# conda install mamba
# mamba install sra-tools fastqc trimmomatic multiqc curl 
# mamba install bwa samtools picard gatk4
# mamba install snpeff snpsift
```

---

## 2. Genomic Sequencing Data

### Dataset Information
- The sequencing data can be retrieved from publicly available datasets:
  *The identifiers SRR1972917 and SRR1972918 correspond to sequencing runs from a study on the 2014 Ebola virus outbreak in Sierra Leone. These runs are part of the             BioProject PRJNA257197 and were sequenced using the Illumina HiSeq 2500 platform. The data consists of 101 base pair paired-end reads obtained from patient blood             samples through metagenomic RNA sequencing approaches. The resulting sequences belong to the SL3 clade of the Makona variant of Zaire ebolavirus.*
- Study on the 2014 Ebola virus outbreak in Sierra Leone [(PMC4503805)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4503805/)
- **SRA ID:** [SRR1972917](https://www-ncbi-nlm-nih-gov.eres.library.manoa.hawaii.edu/sra/?term=SRR1972917)
- **SRA ID:** [SRR1972918](https://www-ncbi-nlm-nih-gov.eres.library.manoa.hawaii.edu/sra/?term=SRR1972918)
- We’ll only use 100k reads
- Use `fastq-dump` to download paired-end reads:

```
mkdir Germline_Variants
cd Germline_Variants/
```
```
# Download raw sequencing data
fastq-dump --split-files -X 100000 SRR1972917
fastq-dump --split-files -X 100000 SRR1972918
```

---

## 3. QC raw reads

### Run FASTQC and Trimmomatic for Quality Filtering and Trimming 
- Quality control (QC) ensures high-quality sequencing reads by identifying and filtering low-quality sequences.
- Input: Raw FASTQ files.
- Output: Quality-controlled and trimmed FASTQ files.

```
mkdir qc
fastqc *.fastq -o qc/
```

```
# Download adapter file and trim sequences
curl -OL https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa > adapters.fa
trimmomatic PE SRR1972917_1.fastq SRR1972917_2.fastq \
    SRR1972917_trimmed_1.fastq SRR1972917_unpaired_1.fastq \
    SRR1972917_trimmed_2.fastq SRR1972917_unpaired_2.fastq \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```

```
trimmomatic PE SRR1972918_1.fastq SRR1972918_2.fastq \
    SRR1972918_trimmed_1.fastq SRR1972918_unpaired_1.fastq \
    SRR1972918_trimmed_2.fastq SRR1972918_unpaired_2.fastq \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```

```
mkdir qc_trimmed
fastqc *trimmed_*.fastq -o qc_trimmed/
```


---

## 4. Align reads

### Run bwa mem
- Reads are aligned to the reference genome using BWA MEM, which is optimized for mapping short reads.
- Input: Trimmed FASTQ files.
- Output: SAM file with aligned reads.
  
```
bwa mem -R '@RG\tID:SRR1972917\tSM:SRR1972917\tPL:ILLUMINA\tLB:SRR1972917' \
    /home/bqhs/ebola/AF086833.fa SRR1972917_trimmed_1.fastq SRR1972917_trimmed_2.fastq > SRR1972917_raw.sam
```
```
bwa mem -R '@RG\tID:SRR1972918\tSM:SRR1972918\tPL:ILLUMINA\tLB:SRR1972918' \
    /home/bqhs/ebola/AF086833.fa SRR1972918_trimmed_1.fastq SRR1972918_trimmed_2.fastq > SRR1972918_raw.sam
```

---

## 5. Sort alignment

### Run samtools sort
- Sorting the aligned reads is necessary to optimize downstream processing.
- Input: SAM file.
- Output: Sorted BAM file.

```
samtools sort SRR1972917_raw.sam > SRR1972917_sort.bam
samtools sort SRR1972918_raw.sam > SRR1972918_sort.bam
```

---

## 6. Mark duplicates

### Run picard MarkDuplicates
- Input: Sorted BAM file.
- Output: BAM file with marked duplicates.
- Marking duplicates removes duplicate reads originating from PCR amplification, reducing bias.
  
```
picard MarkDuplicates -Xmx50g I=SRR1972917_sort.bam O=SRR1972917_dedup.bam M=SRR1972917_dedup.txt
picard CollectAlignmentSummaryMetrics -Xmx50g INPUT=SRR1972917_dedup.bam OUTPUT=SRR1972917_aln_metrics.txt REFERENCE_SEQUENCE=/home/bqhs/ebola/AF086833.fa
samtools flagstat SRR1972917_dedup.bam
```
```
picard MarkDuplicates -Xmx50g I=SRR1972918_sort.bam O=SRR1972918_dedup.bam M=SRR1972918_dedup.txt
picard CollectAlignmentSummaryMetrics -Xmx50g INPUT=SRR1972918_dedup.bam OUTPUT=SRR1972918_aln_metrics.txt REFERENCE_SEQUENCE=/home/bqhs/ebola/AF086833.fa
samtools flagstat SRR1972918_dedup.bam
```

---


## 7. Index BAM file

### Run samtools index
- Input: Deduplicated BAM file.
- Output: Indexed BAM file.
- Indexing enables fast access to specific genomic regions in BAM files.

```
samtools index SRR1972917_dedup.bam
samtools index SRR1972918_dedup.bam
```

---

## 8. Variant Calling

### Run GATK HaplotypeCaller
- HaplotypeCaller is a tool from the GATK (Genome Analysis Toolkit) used for variant discovery (i.e., finding SNPs and small indels).
- Input: Indexed BAM file.
- Output: GVCF files with called variants.
- HaplotypeCaller identifies potential germline variants in the sequencing data.
  
```
gatk HaplotypeCaller -R /home/bqhs/ebola/AF086833.fa -I SRR1972917_dedup.bam -O SRR1972917.g.vcf -ERC GVCF
gatk HaplotypeCaller -R /home/bqhs/ebola/AF086833.fa -I SRR1972918_dedup.bam -O SRR1972918.g.vcf -ERC GVCF
```

### Why use `-ERC GVCF`?
+ This instructs HaplotypeCaller to output an intermediate GVCF (Genomic VCF) file.
+ Unlike a standard VCF which contains only sites determined to be variant, a GVCF contains information for all sites in the genome (both variant and non-variant), accompanied by reference confidence information.
+ GVCF output is especially useful to perform joint genotyping on multiple samples later.

---

## 9. Combine Variants

### Run GATK CombineGVCFs
- This command merges individual GVCF files into one combined GVCF file. 
- Input: Multiple GVCF files.
- Output: Combined GVCF file.
- Combining multiple variant calls into a single file for joint genotyping.

```
gatk CombineGVCFs -R /home/bqhs/ebola/AF086833.fa -V SRR1972917.g.vcf -V SRR1972918.g.vcf -O combined.g.vcf
```
---

## 10. Genotype Variants

### Run GATK GenotypeGVCFs
- This command uses the combined GVCF file (`combined.g.vcf`) to perform the final joint genotyping step and produce a standard VCF (`combined.vcf`) containing the called variants for all samples. 
- Input: Combined GVCF file.
- Output: Final VCF file with genotyped variants.

```
gatk GenotypeGVCFs -R /home/bqhs/ebola/AF086833.fa -V combined.g.vcf -O combined.vcf
```
### What's the difference between GVCF and VCF?
+ GVCF (Genomic VCF) is an intermediate format containing both variant sites and non-variant sites with reference confidence intervals useful for multi-sample workflows.
+ VCF (Variant Call Format) is the final output that contains only the detected variants.

---


## 11. Variant Filtering

### Run GATK VariantFiltration
### 1st Filtering Step (Site-Level Filtering)
- This command applies variant filtration to the combined.vcf file, marking variants that do not meet quality and depth thresholds.
- The output vcf file contains all variants from  input vcf file, but those that do not meet the filtering criteria are flagged with "lowQualDp".
- You can later remove or ignore these flagged variants depending on your analysis needs.
- Arguments:
  + `-R /home/bqhs/ebola/AF086833.fa`: Specifies the reference genome file (AF086833.fa).
  + `-V combined.vcf`: Input VCF file (combined.vcf).
  + `-O combined.filter1.vcf`: Output VCF file (combined.filter1.vcf) after applying the variant filters.
  + `-filter *QUAL < 30.0 || DP < 10*`  Applies a filter to flag variants with:
    + QUAL < 30.0: Variants with a quality score lower than 30.
    + DP < 10: Variants with a depth (number of supporting reads) lower than 10.
  + `--filter-name lowQualDp`: Labels filtered variants as "lowQualDp" in the VCF file.

```
gatk VariantFiltration -R /home/bqhs/ebola/AF086833.fa -V combined.vcf -O combined.filter1.vcf \
    -filter "QUAL < 30.0 || DP < 10" --filter-name lowQualDp
```

### 2nd Filtering Step (Genotype-Level Filtering)
- This step flags low-confidence genotype calls at the sample level, ensuring that only reliable genotype calls are retained.
- Arguments:
  + `-R /home/bqhs/ebola/AF086833.fa`: Specifies the reference genome file (AF086833.fa).
  + `-V combined.filter1.vcf`: Takes the previously filtered VCF file (combined.filter1.vcf) as input.
  + `-O combined.filter2.vcf`: Outputs a new VCF file (combined.filter2.vcf) after applying genotype filters.
  + `-G-filter "GQ < 20.0"`: Applies a genotype-level filter to flag individual genotype calls where:
      + GQ < 20.0: Genotype quality score is below 20.
  + `-G-filter-name lowGQ`: Labels filtered genotype calls as "lowGQ" in the VCF file.

```
gatk VariantFiltration -R /home/bqhs/ebola/AF086833.fa -V combined.filter1.vcf -O combined.filter2.vcf \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ
```

### Why Perform Both Steps?
+ The first step ensures only high-confidence variant sites are considered.
+ The second step ensures only high-confidence genotypes are retained for each individual sample.
+ This two-step approach is critical in multi-sample studies, where some samples may have poor genotype confidence even if the site passes the variant-level filters.

---

## 12. Genotype Concordance
### Run GATK GenotypeConcordance
- Evaluates the accuracy of genotype calls by comparing them to known truth datasets.
### When should you run Genotype Concordance?
- If you have a known truth dataset (e.g., HapMap, GIAB) and want to evaluate accuracy.
- If you're comparing different variant calling methods or different filtering strategies.
- If you need to ensure consistency in genotyping across samples or sequencing batches
  + Useful in clinical and forensic applications where accurate genotyping is critical.
### When is Genotype Concordance not needed?
- If you are working with de novo variant calls and don't have a truth set for comparison.
- If your dataset lacks matched ground truth genotypes.

```
gatk GenotypeConcordance -CV combined.filter2.vcf -TV /home/bqhs/ebola/ebola-samples.vcf -O SRR1972917.concordance.grp -CS SRR1972917 -TS SRR1972917
```


---

## 13. Variant Annotation
### Run snpEff ann
- SnpEff annotates variants with functional information, including their impact on genes.

```
# mamba install snpeff snpsift
snpEff ann -v -c /home/vedbar_2025/miniconda3/share/snpeff-5.2-1/snpEff.config -s snpeff.html AF086833 combined.filter2.vcf > combined.ann.vcf
```

---

## 14. Extract Variant Fields
### Run snpSift extractFields
- Input: Annotated VCF file.
- Output: Extracted variant fields in tabular format.
- Extracts relevant fields from annotated variant files for downstream analysis.

```
SnpSift extractFields combined.ann.vcf \
    ID CHROM POS REF ALT QUAL DP FILTER \
    ANN[0].GENE ANN[0].GENEID ANN[0].EFFECT ANN[0].IMPACT \
    ANN[0].BIOTYPE ANN[0].HGVS_C ANN[0].HGVS_P \
    GEN[0].GT GEN[0].GQ GEN[0].FT \
    GEN[1].GT GEN[1].GQ GEN[1].FT > combined.txt
```

---

## 15. Class Exercise

### Create Directory
```
mkdir Class
cd Class  
```

### Run gatk HaplotypeCaller
```
gatk HaplotypeCaller -R /home/bqhs/hg38/genome.fa -I /home/bqhs/dna/mother.bam -O mother.g.vcf -ERC GVCF 
gatk HaplotypeCaller -R /home/bqhs/hg38/genome.fa -I /home/bqhs/dna/father.bam -O father.g.vcf -ERC GVCF 
gatk HaplotypeCaller -R /home/bqhs/hg38/genome.fa -I /home/bqhs/dna/son.bam    -O son.g.vcf -ERC GVCF
```

### Run gatk CombineGVCFs
```
gatk CombineGVCFs -R /home/bqhs/hg38/genome.fa -V mother.g.vcf -V father.g.vcf -V son.g.vcf -O family.g.vcf
```

#*Note: Running gatk on big files takes time. So, you can copy vcf files and proceed to downstream analysis
```
cp /home/bqhs/dna/Results/*vcf* ./
```

### Run gatk GenotypeGVCFs
```
gatk GenotypeGVCFs -R /home/bqhs/hg38/genome.fa -V family.g.vcf -O family.vcf
```

### Run gatk Variant Filtration
```
gatk VariantFiltration -R /home/bqhs/hg38/genome.fa -V family.vcf -O family.filter.vcf -filter "QUAL < 30.0 || DP < 10" --filter-name lowQualDp 
grep "lowQualDp" family.filter.vcf
```
### Run gatk Calculate Genotype Posteriors
```
cp /home/bqhs/dna/trio.ped ./
gatk CalculateGenotypePosteriors -R /home/bqhs/hg38/genome.fa -V family.filter.vcf -ped trio.ped -supporting /home/bqhs/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz  -O family.CGP.vcf
```
### gatk Variant Filtration
```
gatk VariantFiltration -R /home/bqhs/hg38/genome.fa -V family.CGP.vcf -O family.CGP.filter.vcf -G-filter "GQ < 20.0" -G-filter-name lowGQ 
grep "lowGQ" family.CGP.filter.vcf
```
### Run gatk GenotypeConcordance
```
# gatk CollectVariantCallingMetrics –I family.CGP.filter.vcf --DBSNP /home/bqhs/hg38/dbsnp_146.hg38.vcf.gz -O family.CGP.filter.vcf.metrics
# gatk GenotypeConcordance -CV [your callset vcf] -TV [truth set vcf] -O [output name] -CS [sample name in your callset] -TS [sample name in truth set]
```
### Run snpEff ann

```
#conda install snpeff 
#conda install snpsift
#snpEff ann
#SnpSift extractFields
 snpEff  -Xmx10g ann hg38  -v -s snpeff.html family.CGP.filter.vcf > family.ann.vcf
```
### Run SnpSift extractFields
```
SnpSift annotate /home/bqhs/hg38/dbsnp_146.hg38.vcf.gz family.ann.vcf > family.ann2.vcf
```
### Getting a spreadsheet of variants with annotations
```
SnpSift extractFields family.ann2.vcf \
ID CHROM POS REF ALT QUAL DP FILTER \
ANN[0].GENE ANN[0].GENEID ANN[0].EFFECT ANN[0].IMPACT \
ANN[0].BIOTYPE ANN[0].HGVS_C ANN[0].HGVS_P \
GEN[0].GT GEN[0].GQ GEN[0].FT \
GEN[1].GT GEN[1].GQ GEN[1].FT \
GEN[2].GT GEN[2].GQ GEN[2].FT > family.txt
```

### Transfer all results to your local machine via [FileZilla](https://filezilla-project.org/download.php) and view them.

---

## Contributing
#### Contributions to improve this pipeline are welcome!


---
