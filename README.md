# Comparing phasing outputs 🧬

Welcome to this repository, where we focus on the essential task of choosing the best phasing strategy for your genomic data. Our aim is to provide a thorough comparison of read-based and pedigree phasing and determine the most effective approach for a given dataset.

We utilize individuals with pedigree information and perform both types of phasing, comparing the results with mendelian inheritance. The workflow of the phasing process starts with called variants in a VCF and continues as follows: 
1) Filter VCFs
2) Phase with WhatsHap (either pedigree or read-based phasing)
3) Filtering to ensure quality
4) ShapeIt for final phasing (population-based)

This repository provides a step-by-step guide to the phasing process, allowing you to easily replicate our results and make informed decisions about the best phasing strategy for your data. Join us on this exploration of genomic data analysis and the crucial role that phasing plays in the process.

### Step 1: Filtering VCFs
GATK Best-practices technical filtering, removing masked regions and sequencing depth filtering.
* Input: VCF file ($vcf) and reference genome ($ref)
* Script: ```1_FilteringVCFs.sh```
* Output: filtered VCF

### Step 2: Phasing with WhatsHap

* Input: filtered VCF file ($input_vcf), BAM file ($bam), PED file ($ped), reference genome ($ref)
* Script: ```2_WhatsHap.sh```
* Output: phased VCF
