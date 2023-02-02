# Define path
myPATH=/soniasarmiento/Phasing/Switch
# Save the code into a variable
Switch=${myPATH}/SHAPEIT5_switch_static

CHR=Super-Scaffold_14
#VAL = file from where Switch computes Mendelian Inheritance and uses as validation
VAL=UNPHASED_PARENTS_AND_OFFSPRING
      # PED file of the samples
      ped=${myPATH}/PEDFile.ped
#EST = VCF file with Phased offspring
EST=PHASED_OffspringONLY
#Output prefix:
OUTPUT_PREFIX=${myPATH}
# Switch Running line
${Switch} --validation ${VAL}_${CHR}.vcf.gz --estimation ${EST}_${CHR}.vcf.gz --region ${CHR} --pedigree ${ped} --output ${OUTPUT_PREFIX}_${CHR}
