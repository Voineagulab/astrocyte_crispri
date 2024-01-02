# CROPSeq

This repository contains analysis scripts for "CRISPRi Screening of Enhancers in Human Primary Astrocytes Identifies Regulatory Circuitry Disrupted in Alzheimer's Disease", by Green, Sutton et al. (2024)



## Directory descriptions
**Raw_data_processing**
Code to process raw sequencing data from *e.g.*, RNA-seq, TT-seq, ATAC-seq, and scRNA-seq.

**CandidateEnhancerSelection_LibraryDesign**
For the prioritisation of 979 open chromatin regions, and the design of sgRNA to silence them

**CRISPRi_screen**
Processing and analyses of the CRISPRi screen in NHAs

**Nanostring**
Analyses of Nanostring nCounter data, used to validate hits in the CRISPRi screen using an independent method

**Functional_annotation**
The use of published data to annotate gene, enhancers, and variants with functional phenotypic information

**Predictive_models**
Using functional genomic properties of enhancers and genes to predict functional enhancer-gene interactions.


**Manuscript_preparation_plots**
Scripts here generate final figures and tables published in the manuscript. They use output generated in preceding scripts.

## Other
Questions should be directed to Irina Voineagu (i.voineagu@unsw.edu.au)
