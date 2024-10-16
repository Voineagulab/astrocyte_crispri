#############
##
## Transcription factor related manuscript plots
##
## @author: Sam Bagot
## @date: 23-11-01
## edited by GJS, to change output directory of plots 
############

setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_TranscriptionFactors/")
library(ggplot2)
library(cowplot)
library("stringr")
library("ggsignif")
library(ggbeeswarm)
library("tidyr")
source("../FinalFigureFunctions.R")

#axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),

#########Functions to write to disk in a traceable way
  
  pdf_tf <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script TF - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_tf <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script TF - ", title, ".txt"))
    print(toPrint)
    sink()
  }


#########Plotting functions
footprintingSignalPlot <- function (df, title = "") {
  p1 <- ggplot(df, aes(x = pos, y = values, colour = regions, group = regions)) +
    geom_line(alpha = 0.8) +
    standard_theme + theme(legend.margin=margin(-10), legend.title = element_blank()) +
    scale_colour_manual(values = c(rev(pals$Hits_Darker), pals$Primary_Darker[4]), guide = guide_legend(ncol = 1)) +
    #scale_fill_manual(values = c(NSColour$Fill, HitColour$Fill)) +
    scale_x_continuous(expand = c(0,0)) + 
    labs(y = "ATAC Signal", x = "Distance from\nmotif centre (bp)", title = title)
  return(p1)
}

getMaxTFs <- function (test_boundMatrix, x, sig = T ) {
  bound_sig <- test_boundMatrix
  if (sig == T) {
    bound_sig <- test_boundMatrix[test_boundMatrix$LOG10FDR  > -log10(0.05),]
  }
  maxTFs_P <- setNames(aggregate(by = list(bound_sig[,x]), bound_sig[,"LOG10FDR"], max), c(x, "LOG10FDR"))
  max_TFs_P <- merge(maxTFs_P, bound_sig[,c(x, "LOG10FDR", "concat_names")])
  maxTFs_OR <- setNames(aggregate(by = list(bound_sig[,x]), bound_sig[,"Odds Ratio"], max), c(x, "Odds Ratio"))
  max_TFs_OR<- merge(maxTFs_OR, bound_sig[,c(x, "Odds Ratio", "concat_names")])
  maxTFs <- merge(max_TFs_P, max_TFs_OR, by = x)
  colnames(maxTFs)[colnames(maxTFs) %in% c("concat_names.x","concat_names.y")] <- c("concat_names_P", "concat_names_OR")
  maxTFs$TopORandTopP <- paste0(maxTFs$concat_names_P,"/",maxTFs$concat_names_OR)
  return(maxTFs)
}

# familyTFPlots <- function (test_boundMatrix, x = "JASPAR_Family", y = "LOG10FDR", size = "",expand_x = c(0.1, 0.08)) {
#   #Remove OR if below pvalue
#   maxTFs <- getMaxTFs(test_boundMatrix, x)
#   if (y != "LOG10FDR") {
#     p1 <- geom_hline(yintercept=1, colour="grey" )
#     y.title <- y
#     text.label <- "concat_names_OR"
#   } else {
#     p1 <- geom_hline(yintercept=-log10(0.05), colour="grey", linetype = "dashed" )
#     y.title <- "-Log10FDR"
#     text.label <- "concat_names_P"
#   }
#   p1 <- ggplot(test_boundMatrix, aes(x = .data[[x]], y=  .data[[y]], size = .data[[size]],colour = .data[[x]])) + 
#     geom_point(alpha = 0.7) + #, colour = pals$One
#     standard_theme + theme(panel.border = element_blank(), legend.box.margin=margin(-10,-10,-10,-10),
#                            axis.text.x = element_text(vjust = 0.5, hjust=1,  angle = 90)) + #geom_beeswarm()
#     geom_text(data=maxTFs,fontface = "bold",  size = 2.5,nudge_y = 1.25, angle = 0,  #nudge_x = 0.25, 
#               aes(x = .data[[x]], y=  .data[[y]], label=.data[[text.label]])) + 
#     scale_size_continuous(range = c(0.5,4)) +
#     scale_x_discrete(expand = expansion(add = c(0.7, 0.7)), limits=rev) + scale_y_continuous(expand = expansion(mult = c(0.00, 0.2))) + coord_flip() + guides(colour = "none") +
#     labs(x = sub("_"," ",x), y = y.title) + p1
#   return(p1)
# } 

# familyTFPlots_GJS(test_boundMatrix, x = "JASPAR_Class", y = "LOG10FDR", size = "Odds Ratio", expand_x = c(0.05,0.05))

familyTFPlots_GJS <- function (test_boundMatrix, x = "JASPAR_Family", y = "LOG10FDR", size = "",expand_x = c(0.1, 0.08)) {
  #Remove OR if below pvalue
  maxTFs <- getMaxTFs(test_boundMatrix, x)
  if (y != "LOG10FDR") {
    p1 <- geom_hline(yintercept=1, colour="grey" )
    y.title <- y
    text.label <- "concat_names_OR"
  } else {
    p1 <- geom_hline(yintercept=-log10(0.05), colour="grey", linetype = "dashed" )
    y.title <- "-Log10 FDR"
    text.label <- "concat_names_P"
  }
  p1 <- ggplot(test_boundMatrix, aes(x = .data[[x]], y=  .data[[y]], size = .data[[size]],colour = .data[[x]])) + 
    geom_point(alpha = 0.7) + #, colour = pals$One
    standard_theme + theme(panel.border = element_blank(), legend.box.margin=margin(-10,-10,-10,-10),
                           legend.position = c(0.8, 0.3), # axis.title.y = invis, 
                           axis.text.y = element_text(size = 6)) + #geom_beeswarm()
    geom_text(data=maxTFs,fontface = "bold",  size = 2.5,nudge_y = 1.25, angle = 0,  #nudge_x = 0.25, 
              aes(x = .data[[x]], y=  .data[[y]], label=.data[[text.label]])) + 
    guides(size = guide_legend(nrow = 3)) +
    scale_size_continuous(range = c(0.5,4)) +
    scale_x_discrete(expand = expansion(add = c(0.7, 0.7)), limits=rev) + scale_y_continuous(expand = expansion(mult = c(0.00, 0.2))) + coord_flip() + guides(colour = "none") +
    labs(x = sub("_"," ",x), y = y.title) + p1
  
  p1
  # return(p1)
} 

plotTFcounts <- function(TF_counts, log2 = F, y.title = "Number of TF binding sites") {
  if (log2 == T) {
    offset <- 1
    TF_counts$`TF counts` <- TF_counts$`TF counts` + offset
  }
  plot <- ggplot(TF_counts, aes( x = category, y = .data[["TF counts"]], colour = category)) + 
    geom_quasirandom() + standard_theme + #geom_quasirandom() +
    theme(panel.border = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) +
    scale_colour_manual(values = rev(pals$Hits), guide = guide_legend(ncol = 1))  +
    geom_signif(comparisons = list(c("Inactive candidates", "Functional enhancers")), vjust = -0.05, colour = pals$One, textsize = 2) +#
    stat_summary(colour = rep(rev(pals$Hits_Darker), length(unique(TF_counts$name))), fun = mean) +
    labs(y = y.title) #
  if (log2 == T) {
    plot <- plot + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)+ offset)
  } else {
    plot <- plot + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  }
  return(plot)
}

#######
#Footprinting data
res <- readRDS(file = "/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/Aggregated_Signal.rds")
res <- lapply(res, function (x) {
  x <- as.data.frame(x)
  x[,"regions"] <-  as.character(x[,"regions"])
  x[x[,"regions"] == "Hits", "regions"] <- "Functional enhancers"
  x[x[,"regions"] == "Non-Hits", "regions"] <- "Inactive candidates"
  x[x[,"TFBS"] == "Bound_TFs", "TFBS"] <- "Bound TFs"
  x[x[,"TFBS"] == "Unbound_TFs", "TFBS"] <- "Unbound TFs"
  return(x)
})
bound_unbound<- rbind(res$Aggregate_Hits_vs_Non_Hits_Bound_Only, res$Aggregate_Hits_vs_Non_Hits_Unbound_Only)

#######
#load TF bound summaries
TF_counts <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/TOBIAS_TF_counts.csv")
TF_counts <- TF_counts[,colnames(TF_counts) %in% c("Exp_Bound_TF_counts", "Exp_Unbound_TF_counts","Hit")]
colnames(TF_counts) <- sub("Exp_(.*)_TF_counts","\\1",colnames(TF_counts))
TF_counts$category <- "Inactive candidates"
TF_counts[TF_counts$Hit,]$category <- "Functional enhancers"
TF_counts$category <- factor(TF_counts$category, levels = c( "Functional enhancers","Inactive candidates"))
TF_counts[,"Fraction Bound"] <- TF_counts$Bound / (TF_counts$Bound + TF_counts$Unbound)
TF_counts <- pivot_longer(TF_counts, c("Bound", "Unbound","Fraction Bound"), values_to = "TF counts" )
TF_counts$name <- factor(TF_counts$name, levels = c("Bound", "Unbound","Fraction Bound"))


#TOBIAS plot aggregate is slightly bugged they only collect a 120 width region when they should grab 121
# pdf("Tobias_Fooprinting_Bound_and_Unbound.pdf", height = 4, width = w_margin *(1/4))
pdf_tf(figNo = "5A", title = "footprinting signal", h = 4, w = w_margin *(1/4))
p1 <- footprintingSignalPlot(bound_unbound) + facet_wrap(vars(TFBS), ncol = 1)
p2 <- p1 + theme(legend.position = "none")
le1 <- get_legend(p1)
plot_grid(p2, le1, nrow = 2, rel_heights = c(1, 0.2))
dev.off()

# pdf("Tobias_Exp_BoundTFs_log2.pdf", width = w_margin *(1/4) , height = 6)
pdf_tf(figNo = "5A", title = "binding counts and fractions", h = 6, w = w_margin *(1/4))
p1 <- plotTFcounts(TF_counts[TF_counts$name != "Fraction Bound",], y.title =  "Number of TF binding sites",log2 = T)  + facet_wrap(vars(name), scales = "free_y",ncol =1)
p2 <- plotTFcounts(TF_counts[TF_counts$name == "Fraction Bound",], y.title = "Fraction Bound") + facet_wrap(vars(name), scales = "free_y",ncol =1)
p2 <- p2 + theme(legend.position = "none")
le1 <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
plot_grid(p1, p2, le1,  ncol =1, rel_heights = c(1, 0.5,0.2) ) #, le1, 0.2
dev.off()

#Combined footprinting & Expressed TFs

p1 <- footprintingSignalPlot(bound_unbound) + facet_wrap(vars(TFBS), nrow = 1)
footprintPlot <- p1 + theme(legend.position = "none")
le1 <- get_legend(p1)

bound_counts_plot <- plotTFcounts(TF_counts[TF_counts$name != "Fraction Bound",], y.title =  "Number of TF binding sites",log2 = T) + facet_wrap(vars(name),nrow =1)
bound_frac_plot <- plotTFcounts(TF_counts[TF_counts$name == "Fraction Bound",], y.title = "Fraction Bound") + facet_wrap(vars(name), scales = "free_y",nrow =1)
bound_frac_plot <- bound_frac_plot + theme(legend.position = "none", axis.title.y = element_blank(), axis.ticks.x = element_blank())
bound_counts_plot <- bound_counts_plot + theme(legend.position = "none", axis.ticks.x = element_blank())

# pdf("Tobias_combined_bound_plot.pdf", width = w_margin *(3.5/8) , height = 3)
pdf_tf(figNo = "5A", title = "Combined", h = 3, w = w_margin *(3.5/8))
plot_grid(footprintPlot, le1,bound_counts_plot,bound_frac_plot, nrow = 2, ncol =2, rel_widths = c(1, 0.6), rel_heights = c(1, 0.8))
dev.off()

########
#Tobias Class plots
test_boundMatrix <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/JASPAR_group_AnnotatedTFs.csv")
colnames(test_boundMatrix)[colnames(test_boundMatrix) == "Odds.Ratio"] <- "Odds Ratio"
#Alter names so they fit on the plot
test_boundMatrix$JASPAR_Class <- gsub("Nuclear receptors", "Nuc. recep.", test_boundMatrix$JASPAR_Class) # %>%
  # gsub("binding", "", .)# test_boundMatrix$JASPAR_Class <- gsub("Nuclear receptors", "Nuc. recep.", test_boundMatrix$JASPAR_Class) %>%
  # gsub("binding", "", .)

# pdf("Tobias_ClassFTPlots_P.pdf", width = w_margin *(4.5/8) , height = 4)
pdf_tf(figNo = "SF7A", title = "JASPAR class enrichments", h = 3.5, w = 3.6)
# familyTFPlots(test_boundMatrix, x = "JASPAR_Class", y = "LOG10FDR", size = "Odds Ratio", expand_x = c(0.05,0.05))
familyTFPlots_GJS(test_boundMatrix, x = "JASPAR_Class", y = "LOG10FDR", size = "Odds Ratio", expand_x = c(0.05,0.05))
dev.off()


