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

#Subsetting the selected variants
gatk SelectVariants \
	-R $ref \
	-V $intermediateVCF \
	-O $filteredVCF \
	--exclude-filtered true \
	--exclude-non-variants
	
	
# Remove masked regions
bcftools view -O z -o $outputVCF -R maskedregions_150_90.bed \
$filteredVCF.vcf.gz

#Index bcf
tabix $outputVCF

# Filter variants for sequencing depth (DP)
gatk VariantFiltration \
-V ${outvcf} \
-G-filter "DP > $cutoffhigh" \
-G-filter-name 'GDPhigh' \
-G-filter "DP < $cutofflow" \
-G-filter-name 'GDPlow' \
-O  ${outvcf2} \
--set-filtered-genotype-to-no-call true
