## Check expression on gene strands

run_fc <- function(saf, saf_dir, stranded = 1) {
      
      write.table(saf, saf_dir, sep = "\t", row.names = FALSE, quote = FALSE)
      
      featureCounts(c(bam_tt, bam_rnaseq), 
                    annot.ext = saf_dir, isGTFAnnotationFile = FALSE, # a SAF of filtered bins
                    useMetaFeatures = FALSE, # feature-level, rather than pooling features-to-genes
                    allowMultiOverlap = TRUE, # a change versus Irina
                    isPairedEnd = TRUE, # self-explanatory
                    nthreads = 8, # self-explanatory 
                    strandSpecific = stranded, 
                    minMQS = 10, # minimum quality score at either end of the pair
                    checkFragLength = FALSE, 
                    countMultiMappingReads = FALSE, 
                    requireBothEndsMapped = FALSE,                              
                    countChimericFragments = FALSE)
      
      # the change in this versus Irina is that allowMultOverlap = TRUE. Because some of our features are overlapping, this is useful
      # further, because countMultiMappingReads = FALSE, we know they aren't multimappers
    }


gtf <- "/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/cellranger_refdata-gex-GRCh38-2020-A/genes/genes.gtf"
bam_tt <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10844A1.markdup.sorted.bam"
bam_rnaseq <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10856A1.markdup.sorted.bam"


geneCheck_sense <- featureCounts(c(bam_tt, bam_rnaseq), 
                                    annot.ext = gtf, 
                                    isGTFAnnotationFile = TRUE, 
                                    GTF.attrType = "gene_name",
                                    useMetaFeatures = TRUE, # feature-level, rather than pooling features-to-genes
                                    allowMultiOverlap = TRUE, # a change versus Irina
                                    isPairedEnd = TRUE, # self-explanatory
                                 countReadPairs = TRUE,
                                    nthreads = 8, # self-explanatory 
                                    strandSpecific = 1, 
                                    minMQS = 10, # minimum quality score at either end of the pair
                                    checkFragLength = FALSE, 
                                    countMultiMappingReads = FALSE, 
                                    requireBothEndsMapped = FALSE,                              
                                    countChimericFragments = FALSE)

geneCheck_antisense <- featureCounts(c(bam_tt, bam_rnaseq), 
                                    annot.ext = gtf, 
                                    isGTFAnnotationFile = TRUE, 
                                    GTF.attrType = "gene_name",
                                    useMetaFeatures = TRUE, # feature-level, rather than pooling features-to-genes
                                    allowMultiOverlap = TRUE, # a change versus Irina
                                    isPairedEnd = TRUE, # self-explanatory
                                    countReadPairs = TRUE,
                                    nthreads = 8, # self-explanatory 
                                    strandSpecific = 2, 
                                    minMQS = 10, # minimum quality score at either end of the pair
                                    checkFragLength = FALSE, 
                                    countMultiMappingReads = FALSE, 
                                    requireBothEndsMapped = FALSE,                              
                                    countChimericFragments = FALSE)

x <- list(TTseq = data.frame(Sense = geneCheck_sense$counts[,"GOK10844A1.markdup.sorted.bam"],
                             Antisense = geneCheck_antisense$counts[,"GOK10844A1.markdup.sorted.bam"]),
          RNAseq = data.frame(Sense = geneCheck_sense$counts[,"GOK10856A1.markdup.sorted.bam"],
                              Antisense = geneCheck_antisense$counts[,"GOK10856A1.markdup.sorted.bam"]))
x <- lapply(x, function(y) {
  y <- log2(y+1)
  y <- y[unique(res.final$Gene),]
  y$Density <- get_density(y$Sense, y$Antisense, n = 250)
  return(y)
})
x <- do.call("rbind", x)
x$Assay <- splitter(rownames(x), "\\.", 1)

pdf(file = "../../../FullScale/Results/Scratchspace/Sense versus Antisense TTseq.pdf", height = 3.5, width = 8)
ggplot(x, aes(x = Sense, y = Antisense, colour = Density)) +
  geom_point() +
  theme_bw() +
  theme(panel.border = invis, axis.line.x = element_line()) +
  facet_wrap(~Assay, nrow = 1, strip.position = "left") +
  scale_colour_viridis_c() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Mapping to Sense Strand in GTF", y = "Mapping to Antisense Strand in GTF") 
dev.off()


