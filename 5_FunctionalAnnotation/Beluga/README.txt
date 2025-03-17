1_KnownVariants.R Extracts SNPs located within a 1 kb window of the tested enhancers, prepares the data for disease-impact scores (DIS) prediction using Beluga https://hb.flatironinstitute.org/sei/ , processes the output of the online tool, and carries out statistical tests comparing the DIS scores between hit and non-hit enhancers.

2_CreateVCF.R Creates VCFs for all possible sequence variants at each nucleotide position for the 145 functional enhancers for In Silico Mutagenesis analysis using https://hb.flatironinstitute.org/deepsea/?analysis=insilico

3_BelugaVariants.R Processes the Beluga outputs for the 145 functional enhancers

4.1_BedToBigWig.R & 4.2_BedToBigWig.sh Create BigWig files for visualization in genome browsers

5.ISM_Beluga_CCL2.R Compares Disease impact score (DIS) values between SNPs within 1kb of functional hits and powered non-hit and explore the in-silico mutagenesis data  generated for Enh427 (which regulates the expression of CCL2)

getFasta.sh Extract FASTA sequences for enhancer regions

#The usage for running Beluga disease variants analyses was the following:
First run getFasta.sh to extract FASTA sequences from the HG38 genome. 
Next, execute createVCFs.R to generate VCF files, ensuring positions are adjusted by +1 (posplus1) to match Beluga’s reference genome format. 
Then, upload the VCF files to the Beluga web tool, selecting Model: Beluga and Genome: HG38. BelugaVariants.R to parse the disease impact scores (DIS), merge results with enhancer data, and perform downstream analyses. 


