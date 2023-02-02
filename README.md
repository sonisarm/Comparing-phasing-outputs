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
Filtering for minor allele cound (MAC), missing data and excluding indels.
* Input: phased VCF ($phased_vcf)
* Script: ```3_Post-WhatsHap_filtering.sh```
* Output: filtered vcf ($filter2.recode.vcf.gz)

### Step 3: ShapeIt phasing
In this step we phase our genetic data with [ShapeIt4](https://www.nature.com/articles/s41467-019-13225-y) software already phased using WhatsHap because it improves accuracy of haplotype inference and resolves ambiguities in the data, resulting in a more accurate representation of haplotype information for downstream applications. Specifically, ShapeIt uses information from multiple individuals to perform [population-based phasing](https://academic.oup.com/bioinformatics/article/35/14/i242/5529122). It is important to be aware of its limitations when dealing with related individuals, as it may result in a biased outcome. To avoid this, it is recommended to phase related samples with a reference panel of known unrelated individuals.



ShapeIt software enables you to refine genetic data that has already been phased using WhatsHap or another phasing tool. The refinement process performed by ShapeIt involves using the information from multiple individuals to improve the accuracy of haplotype inference and resolve any remaining ambiguities in the phased data. The result is a more accurate and robust representation of haplotype information, which can be used for a variety of downstream applications, such as imputation, association studies, and haplotype-based analyses. By incorporating ShapeIt into your phasing workflow, you can ensure that your results are based on the most accurate and reliable data possible.

**IMPORTANT**: Shapeit is a powerful software for haplotype phasing, however it is important to be aware of its limitations when dealing with related individuals. As the software uses population-based phasing that takes into account the genotypes of multiple individuals, it may result in a biased outcome if there are related individuals present in the sample. To avoid this issue, it is recommended to phase related samples with a reference panel that includes individuals that have been previously phased with Whtashap and Shapeit and are known to be unrelated. To identify unrelated samples in your dataset, you can follow the guide ```FindingUnrelated.md```.

* Input: filtered VCF ($)
* Script: ```4_ShapeIt.sh```
* Output: phased vcf ()





