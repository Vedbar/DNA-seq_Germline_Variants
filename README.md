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
7. [Index BAM file](#7-index-bam-fle)
8. [Variant Calling](#8-variant-calling)
9. [Combine Variants](#9-combine-variants)
10. [Genotype Variants](#10-genotype-variants)
11. [Variant Filtering](#11-variant-filtering)
12. [Genotype Concordance](#12-genotype-concordance)
13. [Variant Annotation](#13-variant-annotation)
14. [Extract Variant Fields](#14-extract-variant-fields)

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
    trimmed_1.fastq unpaired_1.fastq \
    trimmed_2.fastq unpaired_2.fastq \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```
```
mkdir qc_trimmed
fastqc trimmed_*.fastq -o qc_trimmed/
```


---

## 4. Align reads

### Run bwa mem
- Reads are aligned to the reference genome using BWA MEM, which is optimized for mapping short reads.
- Input: Trimmed FASTQ files.
- Output: SAM file with aligned reads.
  
```
bwa mem -R '@RG\tID:SRR1972917\tSM:SRR1972917\tPL:ILLUMINA\tLB:SRR1972917' \
    /home/bqhs/ebola/AF086833.fa trimmed_1.fastq trimmed_2.fastq > SRR1972917_raw.sam
```

---

## 5. Sort alignment

### Run samtools sort
- Sorting the aligned reads is necessary to optimize downstream processing.
- Input: SAM file.
- Output: Sorted BAM file.

```
samtools sort SRR1972917_raw.sam > SRR1972917_sort.bam
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

---

## 7. Index BAM file
### Run samtools index
- Input: Deduplicated BAM file.
- Output: Indexed BAM file.
- Indexing enables fast access to specific genomic regions in BAM files.

```
samtools index SRR1972917_dedup.bam
```

---

## 8. Variant Calling

### Run GATK HaplotypeCaller
- Input: Indexed BAM file.
- Output: GVCF files with called variants.
- HaplotypeCaller identifies potential germline variants in the sequencing data.
  
```
gatk HaplotypeCaller -R /home/bqhs/ebola/AF086833.fa -I SRR1972917_dedup.bam -O SRR1972917.g.vcf -ERC GVCF
gatk HaplotypeCaller -R /home/bqhs/ebola/AF086833.fa -I SRR1972918_dedup.bam -O SRR1972918.g.vcf -ERC GVCF
```

---

## 9. Combine Variants

### RUn GATK CombineGVCFs
- Input: Multiple GVCF files.
- Output: Combined GVCF file.
- Combining multiple variant calls into a single file for joint genotyping.

```
gatk CombineGVCFs -R /home/bqhs/ebola/AF086833.fa -V SRR1972917.g.vcf -V SRR1972918.g.vcf -O combined.g.vcf
```
---

## 10. Genotype Variants

### Run GATK GenotypeGVCFs
- Input: Combined GVCF file.
- Output: Final VCF file with genotyped variants.
- This step determines the most probable genotypes for each variant.

```
gatk GenotypeGVCFs -R /home/bqhs/ebola/AF086833.fa -V combined.g.vcf -O combined.vcf
```

---


## 11. Variant Filtering

### Run GATK VariantFiltration
- Filtering removes low-confidence variants based on quality scores and read depth.
- 
```
gatk VariantFiltration -R /home/bqhs/ebola/AF086833.fa -V combined.vcf -O combined.filter1.vcf \
    -filter "QUAL < 30.0 || DP < 10" --filter-name lowQualDp

gatk VariantFiltration -R /home/bqhs/ebola/AF086833.fa -V combined.filter1.vcf -O combined.filter2.vcf \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ
```

---

## 12. Genotype Concordance
### Run GATK GenotypeConcordance
- Evaluates the accuracy of genotype calls by comparing them to known truth datasets.

```
gatk GenotypeConcordance -CV combined.filter2.vcf -TV /home/bqhs/ebola/ebola-samples.vcf -O SRR1972917.concordance.grp -CS SRR1972917 -TS SRR1972917
```


---

## 13. Variant Annotation
### Run snpEff ann
- SnpEff annotates variants with functional information, including their impact on genes.

```
mamba install snpeff snpsift
snpEff ann AF086833 -v -c /home/vedbar/miniconda3/share/snpeff-5.0-0/snpEff.config -s snpeff.html combined.filter2.vcf > combined.ann.vcf
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

## Contributing
#### Contributions to improve this pipeline are welcome!


---
