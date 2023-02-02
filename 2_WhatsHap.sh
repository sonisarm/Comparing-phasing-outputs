#INPUTS:
#input VCF is the filtered VCF with unphased haplotypes. It should contain ONLY the individuals to be phased
#input BAM contains reads of individuals to be phased - must contain at least the individuals that are being phased.

# Read-base phasing
whatshap phase \
     --ped $ped \
     -o $phased.vcf \
     --reference $ref \
     --indels \
     $input.vcf.gz \
     $bam.bam

     # Zip VCF file
     bcftools view  $phased.vcf -O z -o $phased_zipped.vcf.gz
     #Index VCF file
     bcftools index $phased_zipped.vcf.gz

# ADDITIONAL INPUT:
#ped: PED file containing pedigree information of individuals (Offspring ID | Dad ID | Mom ID | Sex). This file is required for pedigree phasing but not for read-base phasing.
whatshap phase \
     -o $phased.vcf \
     --reference $ref \
     --indels \
     $input.vcf.gz \
     $bam.bam

     # Zip VCF file
     bcftools view  $phased.vcf -O z -o $phased_zipped.vcf.gz
     #Index VCF file
     bcftools index $phased_zipped.vcf.gz
