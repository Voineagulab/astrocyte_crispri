


## Setup
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  sample.colours <- carto_pal(7, "Bold")
  names(sample.colours) <- samples

  
## Read in raw cellranger count outpu
 raw <- lapply(samples, function(j) {
   print(j)
   Read10X(data.dir = paste0("../../Data/Sequencing/FinalCount/", j,"/outs/raw_feature_bc_matrix/"))
 })
 
 names(raw) <- samples  
 
## Plot the output of barcodeRanks, to infer the true UMI threshold
 # function
 knee.plot <- function(sample, top.ranks = 1e5, lower.umi = 1000) {
   x <- barcodeRanks(raw[[sample]], lower = lower.umi)
   
   y <- x[which(x$rank < top.ranks),]
   # p <- data.frame(Rank = y$rank, Total = x$total)
   xlab <- paste0("Barcode Rank\n",
                  length(which(y$total > y@metadata$knee)), "/", length(which(y$total > y@metadata$inflection)), " cells > Knee/Inflection") 
   
   qplot(x = y$rank, y = y$total, geom = "line") +
     scale_y_continuous(trans = "log10", labels = comma, breaks = c(100, 1000, 5000, 10000, 25000, 50000, 100000)) +
     scale_x_continuous(trans = "log10", labels = comma, breaks = c(10, 100, 1000, 5000, 10000, 100000), limits = c(1000, 100000)) +
     theme_bw() +
     theme(panel.border = invis, axis.line = element_line(), axis.text.x = element_text(angle = 45, hjust = 1)) +
     annotate("text", x = 2000, y = y@metadata$knee, colour = "dodgerblue", label = paste0("Knee: ", round(y@metadata$knee)), size = 3) +
     annotate("text", x = 2000, y = y@metadata$inflection, colour = "firebrick1", label = paste0("Inflection: ", round(y@metadata$inflection)), size = 3) +
     # annotate("text", x = 2, y = lower.umi, colour = "black", label = paste0("Background: <", lower.umi)) +
     geom_hline(yintercept = c(y@metadata$inflection, y@metadata$knee, lower.umi), colour = c("firebrick1", "dodgerblue", "black"), linetype = c(1,1,2)) +
     labs(x = xlab, y = "nUMIs", title = sample)
   
 }
 
  # apply
  pdf(file = "Barcode Calling - BarcodeRanks Knee (lower=1000).pdf", height = 10, width = 10)
  p <- lapply(names(raw), knee.plot)
  plot_grid(plotlist = p, ncol = 3)
  dev.off()
  
  pdf(file = "Barcode Calling - BarcodeRanks Knee (lower=500).pdf", height = 6, width = 10)
  p <- lapply(names(raw), knee.plot, lower.umi = 500)
  plot_grid(plotlist = p, ncol = 3)
  dev.off()

## Run empty drops
  # run
  ed <- lapply(raw, emptyDrops, lower = 1000) # ~5min per sample
  save(ed, file = "Barcode Calling - Empty Drops.rda")
  
  ed.inflection <- lapply(raw, function(x) {
    inflection <- barcodeRanks(x, lower = 1000)@metadata$inflection
    return(emptyDrops(x, lower = inflection))
  }) # ~5min per sample
  save(ed.inflection, file = "Barcode Calling - Empty Drops (Inflection).rda")
  
  # overlap with cr calls, and cells already called
  p <- list()
  for (j in names(ed.inflection)) {
    print(j)
    x <- as.data.frame(ed.inflection[[j]]@listData)
    x$Barcode <- paste0(j, "_", colnames(raw[[j]]))
    x$CellRanger <- x$Barcode %in% paste0(j, "_", colnames(Read10X(data.dir = paste0("../../Data/Sequencing/FinalCount/", j,"/outs/filtered_feature_bc_matrix/"))))
    x$Used <- x$Barcode %in% colnames(nha)
    
    p[[j]] <- x
  }
  
  # filter to tested barcodes (i.e., those above the "lower" parameter)
  # p <- lapply(p, function(x) x[-which(is.na(x$LogProb)),]) 
    p <- lapply(p, function(x) x[-which(x$Total < 1000),]) 

  
  # annotate top 10k barcodes per library
  p <- lapply(p, function(x) {
    x$Top10000 <- rank(-x$Total) < 10000
    return(x)
  })
  
  
  # combine
  p <- do.call("rbind", p)
  
  # annotate with significant FDRs
  p$Sig10 <- p$FDR < 0.1
  p$Sig05 <- p$FDR < 0.05
  p$Sig01 <- p$FDR < 0.01
  
  # add library information
  p$Library <- substr(rownames(p), 1, 5)
  
  # add information on threshold
  p$HighExp <- p$Total > 10000

  # save
  write.csv(p, file = "Barcode Calling - Empty Drops (Inflection).csv")

## Compare different methods of calling cells!
  meth.col <- pal_lancet()(4)
  names(meth.col) <- c("CellRanger", "Top10000", "EmptyDrops", "HighExp")
  
  # calls per method per library
  q <- melt(p[,c("CellRanger", "Top10000", "Sig01", "HighExp", "Library")], id.vars = "Library")
  q <- q[which(q$value),]
  levels(q$variable)[3] <- "EmptyDrops"
  
  pdf(file = "Barcode Calling - Total Cells Called Across Methods (Inflection).pdf", height = 7, width = 10)
  ggplot(q, aes(x = variable, fill = variable)) +
    geom_bar() +
    facet_wrap(~Library) +
    NoLegend() +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = meth.col) +
    theme_bw() +
    theme(axis.title.y = invis) +
    labs(y = "Number of Called Cells")
  dev.off()
  
  # expression distribution of calls
  q <- melt(p[,c("CellRanger", "Top10000", "Sig01", "HighExp", "Library", "Total")], id.vars = c("Library", "Total"))
  q <- q[which(q$value),]
  levels(q$variable)[3] <- "EmptyDrops"
  
  pdf(file = "Barcode Calling - UMI Distribution in Cells Called Across Methods (Inflection).pdf", height = 7, width = 10)
  ggplot(q, aes(x = variable, fill = variable, y = Total)) +
    geom_violin(scale = "width") +
    facet_wrap(~Library) +
    NoLegend() +
    coord_flip() +
    scale_y_continuous(expand = c(0,0), trans = "log10", labels = scales::comma) +
    scale_fill_manual(values = meth.col) +
    geom_hline(yintercept = 10000, linetype = 2) +
    theme_bw() +
    theme(axis.title.y = invis) +
    labs(y = "Number of Called Cells")
  dev.off()
  
  
  # relationship to the high expression threshold
  q <- melt(p[,c("CellRanger", "Top10000", "Sig01", "HighExp", "Library")], id.vars = c("Library", "HighExp"))
  q <- q[which(q$value | q$HighExp),]
  levels(q$variable)[3] <- "EmptyDrops"

  pdf(file = "Barcode Calling - High Expression Cells Called Across Methods (Inflection).pdf", height = 7, width = 10)
  ggplot(q, aes(x = variable, fill = HighExp)) +
    geom_bar() +
    facet_wrap(~Library) +
    NoLegend() +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    # scale_fill_manual(values = meth.col) +
    theme_bw() +
    theme(axis.title.y = invis) +
    labs(y = "Number of Called Cells")
  dev.off()
  
  q <- q[which(q$HighExp),] # high expression barcodes...
  q <- q[which(!(q$value)),] # ...which are not called as cells
  
  pdf(file = "Barcode Calling - High Expression Cells NOT Called Across Methods.pdf", height = 7, width = 10)
  ggplot(q, aes(x = variable, fill = variable)) +
    geom_bar() +
    facet_wrap(~Library) +
    NoLegend() +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = meth.col) +
    theme_bw() +
    theme(axis.title.y = invis) +
    labs(y = "Number of High Expression Barcodes Not Called As Cells")
  dev.off()
    
    
## Properties of high-expression cells not called
  # list barcodes
  removed <- p[which((!(p$Sig01) | !(p$CellRanger)) & (p$HighExp)),] # a high expression cell, not called by at least one of CR or EmptyDrops
  removed <- data.frame(Barcode = substr(removed$Barcode, 7, 100),
                        Library = removed$Library,
                        UMI = removed$Total,
                        CellRanger = removed$CellRanger,
                        EmptyDrops = removed$Sig01)
  
  # function
  get.percent <- function(ind, x = exp) { sum(x[ind]) / sum(x) }
  
  qc <- list()
  for (j in 1:nrow(removed)) {
    print(j)
    
    # barcode properties
    lib <- removed$Library[j]
    bar <- removed$Barcode[j]
    umi <- removed$UMI[j]
    exp <- raw[[lib]][,bar]
    
    # qc
    qc[[paste0(lib, "_", bar)]] <- data.frame(nUMI = umi,
                                              Mito = get.percent(ind = grep("^MT-", names(exp))),
                                              RiboS = get.percent(ind = grep("^RPS", names(exp))),
                                              RiboL = get.percent(ind = grep("^RPL", names(exp))),
                                              MALAT1 = get.percent(ind = grep("MALAT1", names(exp))))
    
    

  }
  
  qc <- do.call("rbind", qc)
  qc$RemovedBy <- "." 
  qc$RemovedBy[which(removed$CellRanger & !(removed$EmptyDrops))] <- "EmptyDrops Only"
  qc$RemovedBy[which(removed$EmptyDrops & !(removed$CellRanger))] <- "CellRanger Only"
  qc$RemovedBy[which(!(removed$CellRanger) & !(removed$EmptyDrops))] <- "Both"
  qc$RemovedBy <- factor(qc$RemovedBy)
  levels(qc$RemovedBy) <- paste0(levels(qc$RemovedBy), " (n=", table(qc$RemovedBy), ")")
  
  q <- melt(qc)
  levels(q$variable)[-1] <- paste0(levels(q$variable)[-1], " Fraction")
  
  pdf(file = "Barcode Calling - High Expression Cells NOT Called Across Methods, QC.pdf", height = 7, width = 7)
  ggplot(q, aes(x = RemovedBy, y = value, colour = RemovedBy)) +
    geom_jitter(width = 0.2) +
    stat_summary(fun = mean, geom = "point", colour = "black", shape = "-", size = 10) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    geom_hline(yintercept = 0.1, linetype = 2) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.15), axis.text.x = invis, axis.title = invis)
  dev.off()  
  
  
  
## Compare the two EmptyDrops runs
  dat <- list()
  for (j in samples) {
    x <- data.frame(Library = j,
                           Barcode = paste0(j, "_", rownames(ed[[j]])),
                           nUMI = ed[[j]]$Total,
                           Lower1000 = ed[[j]]$FDR,
                           LowerInflection = ed.inflection[[j]]$FDR)
    
    y <- (is.na(x$Lower1000)); z <- which(is.na(x$LowerInflection))
    
    dat[[j]] <- x[-which((is.na(x$Lower1000) & is.na(x$LowerInflection))),] # not an NA in both runs of EmptyDrops
  }
  
  dat <- do.call("rbind", dat)
  
  # categorise barcodes by which algorithm they're called
  dat$Lower1000[which(is.na(dat$Lower1000))] <- 1
  dat$LowerInflection[which(is.na(dat$LowerInflection))] <- 1
  dat$Sig <- "Neither"
  dat$Sig[which(dat$Lower1000 < 0.01 & dat$LowerInflection < 0.01)] <- "Both"
  dat$Sig[which(dat$Lower1000 > 0.01 & dat$LowerInflection < 0.01)] <- "Inflection Only"
  dat$Sig[which(dat$Lower1000 < 0.01 & dat$LowerInflection > 0.01)] <- "Lower1000 Only"
  
  dat <- dat[-which(dat$Sig == "Neither"),]
  
  # get technical metadata for these cells
  qc <- list()
  for(j in samples) {
    print(j)
    x <- dat$Barcode[which(dat$Library == j)]
    x <- gsub(j, "", x)
    x <- gsub("_", "", x)
    
    y <- raw[[j]][,x]

    libsize <- colSums(y)
    qc[[j]] <- data.frame(Library = j, 
                          Barcode = x, 
                          nUMI = log10(libsize),
                          Mito = colSums(y[grep("^MT-", rownames(y)),]) / libsize,
                          RiboS = colSums(y[grep("^RPS", rownames(y)),]) / libsize,
                          RiboL = colSums(y[grep("^RPL", rownames(y)),]) / libsize,
                          MALAT1 = y["MALAT1",] / libsize)
    
  }
  
  qc <- do.call("rbind", qc)
  qc$Sig <- dat$Sig
  q <- melt(qc, id.vars = c("Library", "Barcode", "Sig"))
  
  pdf(file = "Barcode Calling - Properties of Cells Called When Setting Lower To Inflection or 1000.pdf", height = 7, width = 7)
  ggplot(q, aes(x = Sig, y = value, fill = Library)) +
    geom_violin(scale = "width", colour = "black", draw_quantiles = 0.5) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    geom_hline(yintercept = 0.1, linetype = 2) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_fill_manual(values = sample.colours) +
    theme_bw() +
    labs(y = "Value", x = "Cell Called Using EmptyDrop Paramater") +
    theme(legend.position = c(0.8, 0.15))
  dev.off()
  
  # x <- ed.inflection
  # x <- lapply(x, function(y) {
  #   y$Sig <- y$FDR < 0.01
  #   return(y)
  # })  
  # 
  # lapply(x, function(y) table(y$Sig))  
  