# FILTERING CODE:
#bcftools index $phased.vcf #index in case VCF is not phased

# Remove sites with minor allele count less than 5
vcftools --gzvcf $phased.vcf --mac 5 --recode --recode-INFO-all --out $filter1
bcftools view $filter1.recode.vcf -O z -o $filter1.recode.vcf.gz
# Remove indels and sites that have at least 10% missing genotype data
vcftools --vcf ${filter1}.recode.vcf --remove-indels --max-missing 0.9 --recode --recode-INFO-all --out ${filter2}
bcftools view ${filter2}.recode.vcf -O z -o ${filter2}.recode.vcf.gz
bcftools index ${filter2}.recode.vcf.gz
