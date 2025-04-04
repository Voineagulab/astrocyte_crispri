## Herein I take sgRNAs designed for our 979 candidate enhancers, explore their properties, and select the optimal

################################################################################################################################ #
## Setup ----

rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/")
options(stringsAsFactors = FALSE)


################################################################################################################################ #
## Design guides for long peaks, i.e. those >= 200bp in width ----

## Parameters
s <- 200 # the assumed effective range of repression by dCas9-KRAB around a guide binding site
remove.filtered.guides <- FALSE
distance.between.guides <- 50 # the minimum distance between any two pairs of guides

## Functions
  expand.bp <- function(bp, by = s, end = max(l)) {
    range <- (bp - by) : (bp + by)
    range <- range[which(range > 0 & range <= end)] 
    return(range)
  }
  
  collapse.bp <- function(list) {
    vector <- do.call("c", list)
    vector <- unique(vector)
    vector <- vector[order(vector)]
    return(vector)
  }

## Guides
  long.guides <- read.csv("GuideDesign/sgRNAs_long_all.csv")
  # long.guides <- long.guides[c("sgRNA.Sequence", "Strand.of.sgRNA"),]
  
  # filter out lines with erroneous output
  g <- grep("chr", long.guides$Input) # 5 lines don't pass this criterion, they are mostly junk
  long.guides <- long.guides[g,]
  
  final.guides.long <- list()
  regions <- names(table(long.guides$Input))
  
 ## Now run loop
  for (j in regions) { # test on regions[120]
    print(j)
    
    ## Extract information for the target region
      # subset
      x <- long.guides[which(long.guides$Input == j),]
      
      # some guides were excluded from the ranking system due to, e.g. too close to earlier pick, or 
      
      # get length of target region
      l <- x$Target.Total.Length
      l <- 1:l[1] # the range of basepairs to cover
      
      # get positions of guide binding, and extend by +-sbp
      pos <- x$sgRNA.Cut.Position..1.based.
      
    ## Pick the top n guides and check if the entire region is covered
      n <- 2 # minimum number to check
      fully.tiled <- FALSE
    
      # start by checking the first guide...
      use.these <- pos[1:n] 
      y <- x[1:n,] # y is the subset of x which has guides you want to use
      
      covered <- lapply(use.these, expand.bp)
      covered <- collapse.bp(covered)
      uncovered <- l[which(!(l %in% covered)) ]
      
    ## Now check for additional guides, noting that at least one guide is added
      while(!(fully.tiled)) {
        # look at the n+1 guide
        n <- n + 1
        new <- pos[n]
        
        # check if too close to already chosen guides
        dist <- abs(use.these - new)
        dist <- dist < distance.between.guides
        
        if (any(dist)) {
          # note: cannot use next here, as that can cause problems with meeting while conditions when n == length(pos)
          new.coverage <- 0
        } else { # if not too close, will it help?
          new.range <- expand.bp(new)
          new.coverage <- length(which(new.range %in% uncovered))
        }
        

        # and if it helps... update coverage information
        if (new.coverage > 0) {
          use.these <- c(use.these, new)
          y <- rbind(y, x[n,])
          
          
          covered <- lapply(use.these, expand.bp)
          covered <- collapse.bp(covered)
          
          uncovered <- l[which(!(l %in% covered)) ]
        } 
        
        # iterate if there are still uncovered regions or only 1 guide, escape if not!
        if (length(uncovered) == 0 & length(use.these) > 1) fully.tiled <- TRUE
        
        
        ## The below lines of code were useful when initial n == 1, as it ensures that at least 2 guides are picked
        # # iterate if there are still uncovered regions or only 1 guide, escape if not!
        # if (length(uncovered) == 0 & length(use.these) > 1) fully.tiled <- TRUE
        # 
        # ## if the final guide is evaluated and the end conditions are not met
        # if (!(fully.tiled) & n == length(pos)) {
        #   # then choose the highest-ranged guide that is decently spaced
        #   y <- abs(pos - pos[1]) # distance from all guides to the highest-ranked guide (where pos[1] always == use.these[1])
        #   z <- pos[which(y > distance.between.guides )] 
        #   z <- z[1] # the best guide that is not too close to the highest-ranked guid
        #   
        #   use.these <- c(use.these, z)
        #   
        #   fully.tiled <- TRUE
        # }
        # 
      } # end while
      
    ## Get the relevant guides

      final.guides.long[[j]] <- y
  }
  
  final.guides.long <- do.call("rbind", final.guides.long)
  
  write.csv(final.guides.long, file = "GuideDesign/Final Guides, Long Regions (NEW).csv", row.names = FALSE)
 
    
    
################################################################################################################################ #
## On short peaks ----

## Sefi has provided us with a full list of guides, as well as the top 2 for each region
    
## Read in the Broad algorithms top picks
  short.guides.best <- read.csv("GuideDesign/sgRNAs_short_2_best_sgRNAs.csv")
    
    
## Calculate distances between guides
  short.guide.separation <- list()
  regions <- names(table(short.guides.best$Input))
  
  for (j in regions) {
    x <- short.guides.best$sgRNA.Cut.Position..1.based.[which(short.guides.best$Input == j)]
    short.guide.separation[[j]] <- abs(x[1] - x[2])
    
  }
  plot.data <- do.call("rbind", short.guide.separation)
  plot.data <- as.data.frame(plot.data)
  
  pdf(file = "GuideDesign/Short Peaks - Distance Between Sefi's Guides.pdf", height = 3, width = 8)
  ggplot(plot.data, aes(x = V1)) +
    geom_histogram(colour = "black", fill = "darkorange1", breaks = seq(-1, 140, 5)) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(-1,140), expand = c(0,0), breaks = seq(0, 140, 10)) +
    # scale_x_continuous(breaks = function(x) 5*x) +
    geom_vline(xintercept = 0) +
    theme(panel.border = element_blank()) +
    labs(y = "Count in Bin", x = "Distance between guide pair")
  dev.off()
  
## Calculate distances from peak centre
  ## No calculation necessary - just note that 100bp is the middle!
  pdf(file = "GuideDesign/Short Peaks - Start Position of Sefi's Guides.pdf", height = 3, width = 8)
  ggplot(short.guides.best, aes(x = sgRNA.Cut.Position..1.based.)) +
        geom_histogram(colour = "black", fill = "darkorange1", breaks = seq(-1, 200, 5)) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(-1,201), expand = c(0,0), breaks = seq(0, 200, 10)) +
    # scale_x_continuous(breaks = function(x) 5*x) +
    geom_vline(xintercept = c(0,100), linetype = c(1,2)) +
    theme(panel.border = element_blank()) +
    labs(y = "Count in Bin", x = "Start Position of Guide in 200bp Window")
  dev.off()
  
  
## Based on these, the algorithm is biased in its selection of the top2 peaks. 
## Specifically, high-ranking guides will be excluded if their start position is ~ <5bp or >130bp. this is not relevant to our study of these (narrower) peaks, as the suppression window will still cover this
## I will therefore directly select the top2 guides, ignoring this filter.
  
  # read in
  short.guides <- read.csv("GuideDesign/sgRNAs_short_all.csv")
  
  # for each region, pick the two with the best combined rank
  final.guides.short <- short.guides[which(short.guides$Combined.Rank %in% c(1, 2)),]
  
  # save
  write.csv(final.guides.short, file = "GuideDesign/Final Guides, Short Regions.csv", row.names = FALSE)

  
   
  
