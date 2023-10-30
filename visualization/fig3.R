# DEGs by TIMEPOINT
# MJ Olmo-Uceda
################################################################################
## A) Number of DEGs by time. Split UP and DOWN
# B) Number of significant genes per time
str(nDEGs)

ggplot(nDEGs) +
  geom_bar(aes(x = UP,
               y = time),
           stat = "identity",
           fill = up,
           alpha = 0.2) +
  geom_bar(aes(x = up.lfc1,
               y = time),
           stat = "identity",
           fill = up,
           alpha = 0.4) +
  geom_bar(aes(x = up.lfc2,
               y = time),
           stat = "identity",
           fill = up,
           alpha = 0.5) +
  geom_bar(aes(x = up.lfc3,
               y = time),
           stat = "identity",
           fill = up,
           alpha = 0.7) +
  geom_bar(aes(x = -DOWN,
               y = time),
           stat = "identity",
           fill = down,
           alpha = 0.1) +
  geom_bar(aes(x = -DOWN.lfc1,
               y = time),
           stat = "identity",
           fill = down,
           alpha = 0.2) +
  geom_bar(aes(x = -DOWN.lfc2,
               y = time),
           stat = "identity",
           fill = down,
           alpha = 0.4) +
  geom_bar(aes(x = -DOWN.lfc3,
               y = time),
           stat = "identity",
           fill = down,
           alpha = 0.7) +
  # vertical line in 0
  geom_vline(xintercept = 0, color = 1) +

  coord_flip() +
    
  ylab("Time (hpi)") +
  xlab("DEGs (padj <= 0.05)") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank()
  )


################################################################################
# B) upset with the intersections between each contrast?
names(degs.tp.cele) <- c("2 hpi", "4 hpi", "6 hpi", "8 hpi", "12 hpi", "18 hpi", "22 hpi", 
                         "26 hpi", "32 hpi", "38 hpi", "44 hpi")
upset(UpSetR::fromList(degs.tp.cele),
      nsets = length(degs.tp.cele),
      #nintersects = NA,
      sets = rev(names(degs.tp.cele)),
      point.size = 4,
      sets.x.label = "Number DEGs",
      text.scale = 1,
      line.size = 0.4,
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      keep.order = T,
      #expression = "ColName > 3",
      matrix.dot.alpha = 0.2,
      sets.bar.color = "gray50",
      queries = list(list(query = intersects,
                          params = list("6 hpi", 
                                        "8 hpi", "12 hpi", "18 hpi", 
                                        "22 hpi", "26 hpi", "32 hpi", 
                                        "38 hpi", "44 hpi"),
                          color = "#586974",#B83253",
                          active = T),
                     list(query = intersects,
                          params = list("8 hpi", "12 hpi", "18 hpi", 
                                        "22 hpi", "26 hpi", "32 hpi"),
                          color = "#BAB8D4",#C00384",
                          active = T))
      )

names(all.degs)
################################################################################
# C) Volcano plots per time point
#' plotVolcanoFromResuls
#' 
#' First generate a tibble without restrictions, then do the volcano plot with 
#' some anotations
#'
#' @param results DESeqResults object
#' @param padj.cuttof 
#' @param lfc.cutoff 
#'
#' @return
#' @export
#'
plotVolcanoFromResuls <- function(results,
                                  padj.cuttof = 0.05 ,
                                  lfc.cutoff = 1,
                                  color_up = up,
                                  color_down = down,
                                  color_no = "black",
                                  gene = "symbol"
)
{
  # Formatting the data
  results$diffexpressed <- "NO"
  # labeling if the gene is significante
  results$diffexpressed[results$log2FoldChange > lfc.cutoff & 
                          results$pvalue < padj.cuttof] <- "UP"
  results$diffexpressed[results$log2FoldChange < -lfc.cutoff & 
                          results$padj < padj.cuttof] <- "DOWN"
  
  # assigning colors
  colors <- c(color_up, color_no, color_down)
  names(colors) <- c("UP", "NO", "DOWN")
  
  # Creating a column for the labels, only if they are diff expressed will be
  # shown in the final plot
  results$delabel <- NA
  results$delabel[results$diffexpressed != "NO"] <- results$symbol[results$diffexpressed != "NO"]
  
  # VOLCANO PLOT
  volcanoPlot <- ggplot(data = results, 
                        aes(x = log2FoldChange, 
                            y = -log10(padj),
                            col = diffexpressed, 
                            label = delabel)) +
    geom_point(alpha = 0.8) +
    geom_text_repel(aes(fontface = ifelse(diffexpressed != "NO", "italic", "plain")),
                    max.overlaps = 10) +
    scale_color_manual(values=colors) +
    #geom_hline(yintercept=-log10(padj.cuttof), col="red") +
    #geom_vline(xintercept=c(-lfc.cutoff, lfc.cutoff), col="red") +
    theme_minimal()
  
  return(volcanoPlot)
}

## List with the volcano plots to show as grid
volcanos.tp <- list() 
for (res in names(res.tp)){
  volcano <- plotVolcanoFromResuls(res.tp[[res]],
                        lfc.cutoff = 2) +
    labs(title = res,
         x = expression(log[2]("Fold Change")), y = expression(-log[10]("p.adj"))) +
    theme_classic() +
    theme(legend.position = "none")
  volcanos.tp[[res]] <- volcano
}
grid.arrange(grobs = volcanos.tp,
             nrow = 3,
             ncol = 4)

################################################################################
# Number of enriched terms (wormEnrichr) by TP
################################################################################
tp.wormenrichr %>%
  filter(Adjusted.P.value <= 0.05,
         db %in% c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "KEGG_2019")
  ) %>%
  group_by(db, hpi, change) %>%
  summarise(counts = n()) %>%
  ggplot(aes(x = as.factor(hpi),
             y = ifelse(change == "UP", counts, -counts),
             fill = change)) +
  geom_segment(aes(xend = as.factor(hpi), 
                   yend = 0,
                   color = change), 
               size = 0.6) +  # Draw vertical lines
  geom_point(shape = 21, size = 2.5) +  # Draw circular points
  geom_hline(yintercept = 0) +  # Add a horizontal line at y = 0
  facet_wrap(~ db, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("UP" = "#009B9F", "DOWN" = "#c800c9")) +
  scale_color_manual(values = c("UP" = "#009B9F", "DOWN" = "#c800c9")) +
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(colour = "black"))



degs.2waives <- intersect(all.degs$`12 hpi`, unique(c(all.degs$`22 hpi`,
                                      all.degs$`26 hpi`,
                                      all.degs$`32 hpi`)))

degs.2waives.switch <- c()

for (gene in degs.2waives){
  if(sign(betas[gene, "Groupinfection.Time12"]) != sign(betas[gene, "Groupinfection.Time22"]) |
     sign(betas[gene, "Groupinfection.Time12"]) != sign(betas[gene, "Groupinfection.Time26"]) |
     sign(betas[gene, "Groupinfection.Time12"]) != sign(betas[gene, "Groupinfection.Time32"])){
    print(gene)
    degs.2waives.switch <- c(degs.2waives.switch, gene)
  }
}

for (gene in degs.2waives.switch){
  print(plotCounts(dds.lrt, 
                   gene,
                   intgroup = c("Time","Group"), 
                   returnData = TRUE) %>%
        ggplot(.,
               aes(x = as.numeric(as.character(Time)), 
                   y = count, 
                   color = Group, 
                   group = Group)) + 
        geom_point() + 
        ggtitle(cele_trad[match(gene, cele_trad$cele_id), "symbol"]) +
        stat_summary(fun = mean,
                     geom = "line") +
        scale_y_log10() + 
        ylab("Expression level (log10(DESeq2 normalized counts))") +
        xlab("Time point (hpi)") +
        scale_color_manual(values = c("control" = "black",
                                      "infection" = "#7D5BA6")) + #7A5980
        #ggsci::scale_color_jama() +
        theme_classic())
  Sys.sleep(3)
}

# The changes in "cct-3", "C44B7.11" and "T06A1.5" are not so different to 
heatmap.2(betas[intersect(degs.2waives.switch, rownames(betas)),
                grepl("Groupinfection", colnames(betas))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "row",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          labRow = cele_trad[match(intersect(degs.2waives.switch, rownames(betas)), cele_trad$cele_id), 
                             "symbol"],
          cexRow = 0.6,
          adjCol = 0.5,
          srtRow = 20,
          labCol = "",#paste0(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44), " hpi"),
          cexCol = 0.8, 
          srtCol = 20,
          # Add the significance of the log2FC
          cellnote = p.adj.markers[intersect(degs.2waives.switch, rownames(betas)),
          ],
          notecol = "black", 
          notecex = 1,
          
          #adjCol = 0.5,
          #keysize = 0.5,
          key.title = "")


# snap-29, zip-1, "C49C8.8"
plotCounts(dds.lrt, 
           cele_trad[match("zipt-9", cele_trad$symbol), "cele_id"],
           intgroup = c("Time","Group"), 
           returnData = TRUE) %>%
  ggplot(.,
         aes(x = as.numeric(as.character(Time)), 
             y = count, 
             color = Group, 
             group = Group)) + 
  geom_point() + 
  #ggtitle(cele_trad[match(gene, cele_trad$cele_id), "symbol"]) +
  stat_summary(fun = mean,
               geom = "line") +
  scale_y_log10() + 
  ylab("Expression level (log10(DESeq2 normalized counts))") +
  xlab("Time point (hpi)") +
  scale_color_manual(values = c("control" = "black",
                                "infection" = "#7D5BA6")) + #7A5980
  #ggsci::scale_color_jama() +
  theme_classic()


# Genes involved in TGF-beta signaling pathway
heatmap.2(betas.symbol[c(intersect(c("skr-3", "skr-4", "skr-5", "skr-19", "cul-6", "rho-1", "mpk-1", "lin-35"), 
                                 rownames(betas.symbol))),
                grepl("Groupinfection", colnames(betas.symbol))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "row",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          # labRow = cele_trad[match(intersect(degs.2waives.switch, rownames(betas)), cele_trad$cele_id), 
          #                    "symbol"],
          cexRow = 0.6,
          adjCol = 0.5,
          srtRow = 20,
          labCol = "",#paste0(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44), " hpi"),
          cexCol = 0.8, 
          srtCol = 20,
          # Add the significance of the log2FC
          cellnote = p.adj.markers.symbol[intersect(c("skr-3", "skr-4", "skr-5", "skr-19", "cul-6", "rho-1", "mpk-1", "lin-35"), 
                                                    rownames(betas.symbol)),
          ],
          notecol = "black", 
          notecex = 1,
          
          #adjCol = 0.5,
          #keysize = 0.5,
          key.title = "")
