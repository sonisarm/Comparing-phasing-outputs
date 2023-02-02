#Techfilters SNPs selection
gatk VariantFiltration \
	-R $ref \
	-V $VCF \
	-O $tmp \
	--filter-expression "QD < 2.0" --filter-name "QD2" \
	--filter-expression "FS > 60.0" --filter-name "FS60" \
	--filter-expression "SOR > 3.0" --filter-name "SOR3" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS-8" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRS-12.5" \
	--filter-expression "MQ < 40.0" --filter-name "MQ40"

#Subsetting the selected variants
gatk SelectVariants \
	-R $ref \
	-V $tmp \
	-O $filtered.vcf \
	--exclude-filtered true \
	--exclude-non-variants
	
	
# Remove masked regions
bcftools view -O z -o $output.vcf -R maskedregions_150_90.bed $filtered.vcf

#Index bcf
tabix $output.vcf

# Filter variants for sequencing depth (DP)
gatk VariantFiltration \
-V $output.vcf \
-G-filter "DP > $cutoffhigh" \
-G-filter-name 'GDPhigh' \
-G-filter "DP < $cutofflow" \
-G-filter-name 'GDPlow' \
-O  $output_final.vcf \
--set-filtered-genotype-to-no-call true
