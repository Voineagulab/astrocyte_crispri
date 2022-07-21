pdf(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Pos/Test - Effect of MOI.pdf", height = 2.5, width = maxw ) 
for (j in 1:nrow(de.pos)) {
    print(j)

    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]

    p <- data.frame(Exp = nha@assays$VST@data[gene,used],
                    Guide = sceptre.guide.pooled[guide,],
                    MOI = cut(nha$MOI[used], c(-1, 25, 50, 75, 100, 10000)))
    
    levels(p$MOI) <- c("0-25", "26-50", "51-75", "76-100", "101+")

    p$Guide <- factor(p$Guide)
    levels(p$Guide) <- c("Non-targeting", "Targeting")
    tab <- (table(p$Guide))
    levels(p$Guide) <- paste0(levels(p$Guide), " (n=", tab, ")")
    levels(p$Guide)[2] <- paste0(levels(p$Guide)[2], "\nPer Bin:", paste(table(p$Guide, p$MOI)[2,], collapse = ", "))

    pval <- signif(de.pos[j,"Bonf"], 2)
    fc <- signif(exp(de.pos[j,"LogFC"]), 2)


    print(ggplot(p, aes(y = Exp, x = MOI, fill = MOI, shape = Guide)) +
            geom_violin(draw_quantiles = 0.5, colour = "black", adjust = 2, scale = "width") + # note that the adjust argument specifies the extent of smoothing
            geom_jitter(width = 0.2, size = 1) +
            facet_wrap(~Guide) +
            theme_bw() +
            scale_y_continuous() +
            scale_shape_manual(values = c(NA, 1)) +
            labs(y = "NHA Normalised Expression", x = paste0("MOI\nBonferroni=", pval, ", Relative Expression = ", fc)) +
            # scale_fill_manual(values = c("black", "darkorange1")) +
            scale_fill_carto_d(palette = "Magenta") +
            annotate("text", x = 1.5, y = max(p$Exp) * 0.95, label = gene, size = 4) +
            theme(panel.border = invis, axis.line.y = element_line(),
                  legend.position = "none", axis.ticks.x = invis, panel.grid = invis))

}
dev.off()


# Version 2: a scatterplot
pdf(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Pos/Test - Effect of MOI (Scatter).pdf", height = 2.5, width = 3.5) 
for (j in 1:nrow(de.pos)) {
    print(j)

    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]

    p <- data.frame(Exp = nha@assays$VST@data[gene,used],
                    Guide = sceptre.guide.pooled[guide,],
                    MOI = (nha$MOI[used]+1))
    p <- p[order(p$Guide),]
    
    p$Guide <- factor(p$Guide)
    levels(p$Guide) <- c("Non-targeting", "Targeting")
    tab <- (table(p$Guide))
    levels(p$Guide) <- paste0(levels(p$Guide), " (n=", tab, ")")
    levels(p$Guide)[2] <- paste0(levels(p$Guide)[2], "\nPer Bin:", paste(table(p$Guide, p$MOI)[2,], collapse = ", "))

    pval <- signif(de.pos[j,"Bonf"], 2)
    fc <- signif(exp(de.pos[j,"LogFC"]), 2)


    print(ggplot(p, aes(y = Exp, x = MOI, colour = Guide, alpha = Guide)) +
              geom_point(size = 1) +
              scale_alpha_manual(values = c(0.2, 0.5)) +
              theme_bw() +
              geom_smooth(alpha = 1,se = FALSE, linetype = 2) +
              scale_x_continuous(trans = "log2", breaks = c(2, 5, 10, 25, 50, 100, 200, 400)) +
              scale_colour_manual(values = c("black", "firebrick1")) +
              labs(y = paste0("Normalised ", gene, " Expression"), x = paste0("(MOI+1)\nBonferroni=", pval, ", Relative Expression = ", fc), title = paste0(guide, " (n=", tab[2], ")")) +
              theme(panel.border = invis, axis.line = element_line(),
                    legend.position = "none", axis.ticks.x = invis, panel.grid = invis))

}
dev.off()
