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
# mamba install sra-tools fastqc trimmomatic multiqc curl spades quast
mamba install bwa samtools picard gatk4
```

---

## 2. Genomic Sequencing Data

---

## 3. QC raw reads

---

## 5. Align reads

---

## 6. Sort alignment

---

## 7. Mark duplicates

---

## 8. Index BAM file

---

## 9. Variant Calling

---

## 10. Combine Variants

---

## 11. Genotype Variants

---

## 12. Variant Filtering

---

## 13. Genotype Concordance


---

## 14. Variant Annotation

---

## 15. Extract Variant Fields

---
