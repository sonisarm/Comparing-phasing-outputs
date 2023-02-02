shapeit4 --input $post-whatshap_filtered.vcf \
 -T 24 \
 --region $CHR --sequencing \
 --reference $ref \
 --use-PS 0.0001 \
 --pbwt-depth 8 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
 --output $phased_CHR.vcf

bcftools view $phased_CHR.vcf -O z -o $phased_CHR.vcf.gz
bcftools index $phased_CHR.vcf.gz
