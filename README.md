# Comparing phasing outputs 🧬
Author: Sonia Sarmiento

## Introduction
Welcome to this repository, which focuses on the essential task of choosing the best phasing strategy for your genomic data. Our aim is to provide a thorough comparison of read-based and pedigree phasing and determine the most effective approach for a given dataset.

We utilize individuals with pedigree information and perform both types of phasing, comparing the results with Mendelian Inheritance. The workflow of the phasing process starts with called variants in a VCF and continues as follows: 
1) Filter VCFs
2) Phase with WhatsHap (either pedigree or read-based phasing)
3) Filtering to ensure quality
4) ShapeIt for final phasing (population-based)

After this, you have the final phased VCF files for the comparison.
### Aim
This repository provides a step-by-step guide to the phasing process, allowing you to easily replicate our results and make informed decisions about the best phasing strategy for your data. Join us on this exploration of genomic data analysis and the crucial role that phasing plays in the process.

## Workflow

### Step 1: Filtering VCFs
GATK Best-practices technical filtering, removing masked regions and sequencing depth filtering.
* Input: VCF file ($vcf) and reference genome ($ref)
* Script: ```1_FilteringVCFs.sh```
* Output: filtered VCF

### Step 2: Phasing with WhatsHap
Whatshap is a software that allows you to phase genomic data, either using only [read-based](https://www.biorxiv.org/content/10.1101/085050v2.full.pdf) information or using both [read-based and pedigree](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908360/) information if you provide a PED file. With this tool, you can phase both single nucleotide polymorphisms (SNPs) and structural variants (SVs), making it a valuable resource for exploring genetic variation in your data. The process of phasing involves grouping variants that are likely to be inherited together, providing insight into the haplotype structure of your sample.
* Input: filtered VCF file ($input_vcf), BAM file ($bam), PED file ($ped), reference genome ($ref)
* Script: ```2_WhatsHap.sh```
* Output: phased VCF

### Step 3: Filtering VCFs for quality
Filtering for minor allele cound (MAC), missing data and excluding indels.
* Input: phased VCF ($phased_vcf)
* Script: ```3_Post-WhatsHap_filtering.sh```
* Output: filtered vcf ($filter2.recode.vcf.gz)

### Step 4: ShapeIt phasing
In this step we phase our genetic data with [ShapeIt4](https://www.nature.com/articles/s41467-019-13225-y) software already phased using WhatsHap because it improves accuracy of haplotype inference and resolves ambiguities in the data, resulting in a more accurate representation of haplotype information for downstream applications (e.g. imputation). Specifically, ShapeIt uses information from multiple individuals to perform [population-based phasing](https://academic.oup.com/bioinformatics/article/35/14/i242/5529122). It is important to be aware of its limitations when dealing with related individuals, as it may result in a biased outcome. To avoid this, it is recommended to phase related samples with a reference panel of known unrelated individuals, previously phased with WhatsHap and ShapeIt according to this same workflow. To identify unrelated samples in your dataset, you can follow the guide ```FindingUnrelated.md```.

* Input: filtered VCF ($)
* Script: ```4_ShapeIt.sh```
* Output: phased vcf ()


### Step 5 (optional): Merging phased datasets
If you have both unrelated and related samples, you can merge both datasets using the following code. Moreover, before merging you would have to intersect SNPs between the two VCFs.
* Input: phased vcf (unrelated / related) ($unrelated.vcf.gz/related.vcf.gz)
* Script: *below*
* Output: phased merged vcf ($final.vcf.gz)

```bash
# Keep intersect snps of unrelated
bcftools isec -n=2 -w1 -O z -o $unrelated_intersectsnps.vcf.gz unrelated.vcf.gz related.vcf.gz

# Keep intersect snps of related
bcftools isec -n=2 -w1 -O z -o $related_intersectsnps.vcf.gz related.vcf.gz unrelated.vcf.gz

# Merge datasets
bcftools merge -O z -o $final.vcf $unrelated_intersectsnps.vcf.gz $related_intersectsnps.vcf.gz

# If you split your dataset per chromosome/scaffold, you can concatenate the files
bcftools concat --threads 8 --file-list $vcf_ordered.list -O z -o $final.vcf.gz
bcftools index  $final.vcf.gz

# Example of $vcf_ordered.list : 
# ${myPATH}/final_CHR1.vcf.gz
# ${myPATH}/final_CHR2.vcf.gz
# ${myPATH}/final_CHR3.vcf.gz
```



