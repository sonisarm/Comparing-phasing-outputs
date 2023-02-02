     ref='/work/FAC/FBM/DEE/jgoudet/barn_owl/Common/ref_genome_2020/Tyto_reference_Jan2020.fasta' # reference path
        #Extract scaffold files from runLINE
        scaffold=$1
	echo ${scaffold}
  
  
# Phasing UNRELATED individuals
shapeit4 --input ${post-whatshap_filtered.vcf} \
 -T 24 \
 -R ${chr} --sequencing \
 --pbwt-depth 8 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
 --output $phased.vcf
