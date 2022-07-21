# x <- as.matrix(nha@assays$RNA@counts)
# counter <- 0
# x <- apply(x, 2, function(y) {
#   counter <<- counter + 1
#   print(counter)
#   return(y / (sum(y) / 10^6))
# })

library(scuttle)

log2x <- normalizeCounts(as.matrix(nha@assays$RNA@counts)) # per recommendations, use this function

w <- project_cycle_space(x = log2x,
                         species = "human",
                         gname.type = "SYMBOL")

z <- estimate_cycle_position(x = log2x,
                         species = "human",
                         gname.type = "SYMBOL")
# y <- estimate_cycle_position(log2x)
top <- fit_periodic_loess(z, log2x["TOP2A",], plot = TRUE)
y <- estimate_Schwabe_stage(x = log2x,
                         species = "human",
                         gname.type = "SYMBOL")
            
# we also note that users can approximately relate:
  # 0.5pi to be the start of S stage,   
  # pi to be the start of G2M
  # 1.5pi to be the middle of M stage, 
  # 1.75pi-0.25pi to be G1/G0 stage             

## Plot
pdf(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Tricycle.pdf", height = 7, width = 7)
 qplot(w[,1], w[,2]) + theme_bw() + labs(x = "PC1", y = "PC2")
FeaturePlot(nha, features = "Tricycle") + scale_colour_carto_c(palette = "Tropic") + theme_void()
FeaturePlot(nha, features = "Tricycle") + scale_color_gradientn(colours = c("#009B9E", "#F1F1F1", "#C75DAB", "#A16928", "#009B9E"), breaks = c(0, 0.5, 1, 1.5, 2)*pi) + theme_void()
DimPlot(nha, group.by = "TricycleStage", split.by = "TricycleStage", ncol = 2) + scale_colour_lancet(na.translate = TRUE, na.value = "black") + theme_void()
dev.off()


## To do:
  # plot tricycle score on PCA
  # plot Seurat score on PCA
  w <- as.data.frame(w)
  w$Seurat.S <- nha$S.Score
  w$Seurat.G2M <- nha$G2M.Score
  w$Tricycle <- nha$Tricycle
  w$Schwabb <- nha$TricycleStage
  w$Seurat.Stage <- nha$Phase

  base.plot <- ggplot(w, aes(x = PC1, y = PC2))  +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line = element_line(), legend.position = "right") 
    
  pdf(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Tricycle PCA.pdf", height = 4, width = 3.5)
  base.plot + geom_point(mapping = aes(colour = Seurat.S)) +  scale_colour_carto_c(palette = "Geyser")
  base.plot + geom_point(mapping = aes(colour = Seurat.G2M)) +  scale_colour_carto_c(palette = "Geyser")
  base.plot + geom_point(mapping = aes(colour = Tricycle)) +  scale_colour_carto_c(palette = "Tropic")
  base.plot + geom_point(mapping = aes(colour = Schwabb)) +  scale_colour_locuszoom()
  base.plot + geom_point(mapping = aes(colour = Seurat.Stage)) +  scale_colour_locuszoom()
  dev.off()