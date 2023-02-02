# Comparing phasing outputs 🧬

## Introduction
Welcome to this repository, where we focus on the essential task of choosing the best phasing strategy for your genomic data. Our aim is to provide a thorough comparison of read-based and pedigree phasing and determine the most effective approach for a given dataset.

We utilize individuals with pedigree information and perform both types of phasing, comparing the results with Mendelian Inheritance. The workflow of the phasing process starts with called variants in a VCF and continues as follows: 
1) Filter VCFs
2) Phase with WhatsHap (either pedigree or read-based phasing)
3) Filtering to ensure quality
4) ShapeIt for final phasing (population-based)

This repository provides a step-by-step guide to the phasing process, allowing you to easily replicate our results and make informed decisions about the best phasing strategy for your data. Join us on this exploration of genomic data analysis and the crucial role that phasing plays in the process.

## Workflow

### Step 1: Filtering VCFs
GATK Best-practices technical filtering, removing masked regions and sequencing depth filtering.
* Input: VCF file ($vcf) and reference genome ($ref)
* Script: ```1_FilteringVCFs.sh```
* Output: filtered VCF

### Step 2: Phasing with WhatsHap
Whatshap is a software that allows you to phase genomic data, either using only read-based information or using both read-based information and pedigree information if you provide a PED file. With this tool, you can phase both single nucleotide polymorphisms (SNPs) and structural variants (SVs), making it a valuable resource for exploring genetic variation in your data. The process of phasing involves grouping variants that are likely to be inherited together, providing insight into the haplotype structure of your sample.
* Input: filtered VCF file ($input_vcf), BAM file ($bam), PED file ($ped), reference genome ($ref)
* Script: ```2_WhatsHap.sh```
* Output: phased VCF

### Step 3: Filtering VCFs for quality
* Input: phased VCF ($phased_vcf)
* Script: ```3_FilteredVCFs.sh```
* Output: filtered vcf ($filter2.recode.vcf.gz)
