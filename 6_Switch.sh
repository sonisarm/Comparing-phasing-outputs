# Define path
myPATH=/soniasarmiento/Phasing/Switch
# Save the code into a variable
Switch=${myPATH}/SHAPEIT5_switch_static

#Define chromosome/scaffold
CHR=Super-Scaffold_14

#VAL = file from where Switch computes Mendelian Inheritance and uses as validation
VAL=UNPHASED_PARENTS_AND_OFFSPRING
# PED file of the samples
ped=${myPATH}/PEDFile.ped

#EST = VCF file with Phased offspring
EST=PHASED_OffspringONLY
#Output prefix:
OUTPUT_PREFIX=${myPATH}/outputexample

#### SHAPEIT5 switch tool ###
${Switch} --validation ${VAL}_${CHR}.vcf.gz --estimation ${EST}_${CHR}.vcf.gz --region ${CHR} --pedigree ${ped} --output ${OUTPUT_PREFIX}_${CHR}

# this command produce multiple output file.
# You may want to look at the output named $OUTPUT_PREFIX.sample.switch.txt.gz for Switch Error Rate (SER) per sample. This is a 4 columns file, with col1=sample_id and col4=switch error rate.
