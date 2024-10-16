Run order for Beluga disease variants vcf (all positions in Hits)

../getFasta.sh
createVCFs.R
#Run the Beluga web tool https://hb.flatironinstitute.org/sei/ with 
Model = Beluga, Genome = HG38, and the VCFs posplus1 (position plus1)
BelugaVariants.R

(Also using BelugaVariants.R)
