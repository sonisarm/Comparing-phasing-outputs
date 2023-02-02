# Comparing phasing outputs 🧬

Welcome to this repository, where we focus on the essential task of choosing the best phasing strategy for your genomic data. Our aim is to provide a thorough comparison of read-based and pedigree phasing and determine the most effective approach for a given dataset.

We utilize individuals with pedigree information and perform both types of phasing, comparing the results with mendelian inheritance. The workflow of the phasing process starts with called variants in a VCF and continues as follows: 
1) Filter VCFs
2) Phase with WhatsHap (either pedigree or read-based phasing)
3) Filtering to ensure quality
4) ShapeIt for final phasing (population-based)

This repository provides a step-by-step guide to the phasing process, allowing you to easily replicate our results and make informed decisions about the best phasing strategy for your data. Join us on this exploration of genomic data analysis and the crucial role that phasing plays in the process.

## 1) Filtering VCFs


```console
	#Techfilters SNPs selection
	gatk VariantFiltration \
        	-R $ref \
        	-V $VCF \
        	-O $intermediateVCF \
        	--filter-expression "QD < 2.0" --filter-name "QD2" \
        	--filter-expression "FS > 60.0" --filter-name "FS60" \
        	--filter-expression "SOR > 3.0" --filter-name "SOR3" \
        	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS-8" \
        	--filter-expression "MQRankSum < -12.5" --filter-name "MQRS-12.5" \
        	--filter-expression "MQ < 40.0" --filter-name "MQ40"

	#Subsetting
	gatk SelectVariants \
        	-R $ref \
        	-V $intermediateVCF \
        	-O $outputVCF \
        	--exclude-filtered true \
        	--exclude-non-variants

```
