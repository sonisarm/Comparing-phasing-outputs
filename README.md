# Comparing phasing outputs 🧬
Author: Sonia Sarmiento

### This repository provides a step-by-step guide to the phasing process, allowing you to easily replicate our results and make informed decisions about the best phasing strategy for your data. We here use mainly a HPC cluster with bash scripts and R code at the end for plotting.

## Introduction
Welcome to this repository, which focuses on the essential task of choosing the best phasing strategy for your genomic data. Our aim is to provide a thorough comparison of read-based and pedigree phasing and determine the most effective approach for a given dataset.

We utilize individuals with pedigree information and perform both types of phasing, comparing the results with Mendelian Inheritance. The workflow of the phasing process starts with called variants in a VCF and continues as follows: 
1) Filter VCFs
2) Phase with WhatsHap (either pedigree or read-based phasing)
3) Filtering to ensure quality
4) ShapeIt for final phasing (population-based)

After this, you obtain the final phased VCF files for the comparison.

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

* Input: post-whatshap filtered VCF ($post-whatshap_filtered.vcf), reference genome ($ref)
* Script: ```4_ShapeIt.sh```. More information about the code [here](https://odelaneau.github.io/shapeit4/).
* Output: phased vcf per chromosome / scaffold ($phased_CHR.vcf.gz)

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

### Step 6: Comparing phased VCFs with Mendelian Inheritance
SwitchShapeIt 5 is a code for comparing phased VCF genotypes with Mendelian inheritance to assess phasing accuracy. Here, we use it to compare read-base and pedigree phasing. 

**Important**: This procedure can only be performed on offspring with both unphased parents available for comparison.


* Input: phased vcf from offspring (shapeit phased read-base/pedigree), unphased vcf from parents + offsprings (vcfs after ```1_FilteringVCFs.sh```)
* Script: ```6_Switch.sh```
* Output: Switch error rate between Validation (Unphased - Mendelian Inheritance) and Phased VCF. The ouput can be either per sample (```$OUTPUT_PREFIX.sample.switch.txt.gz```) or per SNP (```$OUTPUT_PREFIX.variant.switch.txt.gz```). These are a 4 columns file, with col1=sample_id (for sample.switch.txt.gz only) and col4=switch error rate.

### Step 7: Plotting results
Once you have the output files, you can plot the switch error rate per sample using an R code.
* Input: ```$OUTPUT_PREFIX.sample.switch.txt.gz``` (one per phasing type)
* Script: ```7_Plotting.R```
* Output: ```Results_Sample_Switch.png```


## Results and discussion
The plot suggests that trio phasing is more accurate than read-base phasing after undergoing phasing with shapeit. This conclusion is supported by the study conducted by [Garg et al., 2016](https://academic.oup.com/bioinformatics/article/32/12/i234/2288955), which similarly found that trio phasing was more accurate after undergoing phasing with WhatsHap. A comparison between phasing with one single trio (parents and offspring) and multiple trios (including additional offspring or grandparents) revealed that multiple trios did not lead to improved results and even resulted in greater discrepancies with Mendelian inheritance in some samples. This is contrary to observations from [Blackburn et al., 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7253450/), who determined that phasing real and simulated data from human chromosome X with PULSAR with information of more offspring was more accurate (lower genotype error rate) than with less ofspring (tested 1-7 'children', Figure 3). To conclude, due to our findings and the high computational demands of running WhatsHap with multiple trios, I decided that trio phasing was the most appropriate and accurate method for the data.

Exampe of resulting plot from comparing phasing types with switch tool:
![](https://github.com/sonisarm/Comparing-phasing-outputs/blob/main/Results_Sample_Switch.png)

