################################################################################
# COMPARISON BIOTIC STRESSES.v3 
################################################################################
# Version curated to remain only experiments in N2 WT
# Working with :
#   - pathogens.df --> genes_symbol x Experiments, values: [-1,0,1]
#      
#   - N2.experiments.estreses.df --> experiments x info (86 experiments)
#
# MJ Olmo-Uceda
# 2023/07/19
################################################################################
# Creating DFs
################################################################################
curated.experiments.estreses.df <- readxl::read_xlsx("results/tables/bioticStressesSelected_experiments.xlsx",
                                                     sheet = "curatedSelected")
colnames(curated.experiments.estreses.df)
# Only experiments in N2
N2.experiments.estreses.df <- curated.experiments.estreses.df %>%
  data.frame() %>%
  filter(C..elegans.strain.used == "N2")
rownames(N2.experiments.estreses.df) <- N2.experiments.estreses.df$id
pathogens.df <- pathogens2_wide_df.d[, N2.experiments.estreses.df$id] 

# How many genes are we going to loos
table(rowSums(abs(pathogens.df)) == 0) # 2575
# FALSE  TRUE 
# 14440  2575 
dim(pathogens.df) # 17015    56
# Removing them
pathogens.df <- pathogens.df[rowSums(abs(pathogens.df)) > 0,]
dim(pathogens.df) # 14440    56

### Factorize Category and microbe
N2.experiments.estreses.df$Category <- factor(N2.experiments.estreses.df$Category,
                                              levels = c("Virus", "Gram-", "Gram+", "Fungi", "Microbiota"))
N2.experiments.estreses.df$Microbe.strain.used.for.exposure.treatment <- factor(N2.experiments.estreses.df$Microbe.strain.used.for.exposure.treatment,
                                                                                levels = unique(N2.experiments.estreses.df$Microbe.strain.used.for.exposure.treatment))
################################################################################
# Visualize content of the db
################################################################################
table(N2.experiments.estreses.df$Category)
# Fungi      Gram-      Gram+ Microbiota      Virus 
#     4         21         18          1         12 
################################################################################
# By CATEGORY of pathogen
################################################################################
ggplot(N2.experiments.estreses.df,
       aes(x = Category,
           fill = Category)) +
  geom_bar(#fill = "steelblue", 
    width = 0.5) +
  coord_polar(#theta = "y"
  ) +
  labs(title = "Pathogen Frequency", 
       x = "Category", 
       y = "Experiments") +
  scale_fill_jama() +
  theme_minimal() +
  ylim(c(0,21)) +
  theme(axis.text = element_text(size = 10))

## By PATHOGEN: no microbiota
################################################################################
N2.experiments.estreses.df %>%
  filter(Category != "Microbiota") %>%
ggplot(.,
       aes(x = fct_rev(Microbe.strain.used.for.exposure.treatment),
           fill = Category)) +
  geom_bar(#fill = "steelblue", 
    width = 0.5) +
  #coord_polar() +
  coord_flip() +
  # scale_y_discrete(breaks = seq(0,15)
  #                  #values = seq(0,15)
  #                  ) +
  labs(#title = "Pathogen Frequency", 
    x = "Pathogen", 
    y = "Experiments") +
  # scale_x_discrete(breaks = seq(0, 12, 2),
  #                  labels = seq(0, 12, 2)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Virus" = pal_jama()(5)[1],
                                "Gram-" = pal_jama()(5)[2] , 
                                "Gram+" = pal_jama()(5)[3],
                                "Fungi" = pal_jama()(5)[5])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        legend.text = element_text(size = 5,
                                   colour = "black"),
        axis.text.y = element_text(face = "italic",
                                   colour = "black",
                                   size = 10),
        axis.text.x = element_text(size = 10,
                                   color = "black"))

table(N2.experiments.estreses.df$Category)
################################################################################
# INTERSECTIONS
################################################################################
## OrV responses
OrV.resp.comp <- rownames(pathogens2_wide_df.d[rowSums(abs(pathogens2_wide_df.d %>%
                                                             dplyr::select(rownames(selected.experiments.estreses.df[selected.experiments.estreses.df$Category == "Virus",])))) > 0,])
# Genes that not contribute to the ORV response because we had to remove Sarkies
# from the analysis
setdiff(OrV.resp.comp, 
          rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                      dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Virus",])))) > 0,])
)

OrV.resp.noSarkies <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                    dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Virus",])))) > 0,])
Bth.resp.v3 <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                  dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Microbe.strain.used.for.exposure.treatment %in% c("Bacillus thuringiensis BT247", "Bacillus thuringiensis DB27")  ,])))) > 0,])
## Gram+ responses
gramNeg.resp.v3 <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                                 dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Gram-",])))) > 0,])
## Gram- responses
gramPos.resp.v3 <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                                 dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Gram+",])))) > 0,])
## Fungi
fungi.resp.v3 <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                               dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Fungi",])))) > 0,])
microbiota.resp.v3 <- rownames(pathogens.df[rowSums(abs(pathogens.df %>%
                                                                    dplyr::select(rownames(N2.experiments.estreses.df[N2.experiments.estreses.df$Category == "Microbiota",])))) > 0,])
## List with responses by category
responses.comparison.v3 <- list("OrV" = OrV.resp.noSarkies,
                             "Gram-" = gramNeg.resp.v3,
                             "Gram+" = gramPos.resp.v3,
                             "Fungi" = fungi.resp.v3#,
                             #"Microbiota" = microbiota.resp.v3
                             )
names(responses.comparison.v3)
venn::venn(responses.comparison.v3,
           #zcolor = pal_jama()(5)[c(1:3,5)], #same category colors
           opacity = 0.6,
           #ellipse = T,
           box = F,
           ilcs = 1,
           sncs = 0,
           ggplot = T,
           cex = 0.8
) +
  theme(#text = element_text(size = 10),
    axis.title = element_text(size = 20))

################################################################################
# VISUALIZING GROUP OF GENES
################################################################################
# COMMON GENES
################################################################################
common.v3 <- Reduce(intersect, responses.comparison.v3)
common.v3 %>%
  writeLines()

table(cele_trad[match(orv.specific.v3, cele_trad$symbol), "biotype"])

heatmap.2(as.matrix(pathogens.df[match(Reduce(intersect, responses.comparison.v3),
                                               rownames(pathogens.df)),]),
          #betas[head(order(res.lrt$padj), 40),-c(1,2)],
          dendrogram = "none",
          Colv = F,
          colsep = match(c("by.virus.Orsay..Chen.",
                           "Pseudomonas.sp..vs..OP50",
                           "Bacillus.megaterium",
                           "by.M..humicola..N2...24h."),
                         N2.experiments.estreses.df$experiment),
          sepcolor = "black",
          sepwidth = c(0.2, 0.2),
          #col = c("red", "black", "green"),
          col = c("#c800c9", "white", "#009B9F"),
          #col = "redgreen",
          trace = "none",
          cexRow = 0.7,
          srtRow = 20,
          srtCol = 40,
          #adjCol = 0.5,
          cexCol = 0.5,
          keysize = 0.4,
          ColSideColors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#79AF97FF", "#6A6599FF"
          )[as.numeric(N2.experiments.estreses.df$Category)],
          #key.title = "log2FC"
)

# How are part of the IPR response
intersect(ipr.genes$symbol, 
          Reduce(intersect, responses.comparison.v3)) %>% writeLines()

heatmap.2(as.matrix(pathogens.df[match(intersect(ipr.genes$symbol, 
                                                 Reduce(intersect, responses.comparison.v3)),
                                       rownames(pathogens.df)),]),
          #betas[head(order(res.lrt$padj), 40),-c(1,2)],
          dendrogram = "none",
          Colv = F,
          colsep = match(c("by.virus.Orsay..Chen.",
                           "Pseudomonas.sp..vs..OP50",
                           "Bacillus.megaterium",
                           "by.M..humicola..N2...24h."),
                         N2.experiments.estreses.df$experiment),
          sepcolor = "black",
          sepwidth = c(0.2, 0.2),
          #col = c("red", "black", "green"),
          col = c("#c800c9", "white", "#009B9F"),
          #col = "redgreen",
          trace = "none",
          cexRow = 0.7,
          srtRow = 20,
          srtCol = 40,
          #adjCol = 0.5,
          cexCol = 0.5,
          keysize = 0.4,
          ColSideColors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#79AF97FF", "#6A6599FF"
          )[as.numeric(N2.experiments.estreses.df$Category)],
          #key.title = "log2FC"
)

"eol-1" %in% ipr.genes$symbol

intersect(orv.specific.v3, ipr.genes$symbol) # only 3
# "F57G4.11" "ZC196.3"  "pals-11"

################################################################################
# IPR related
################################################################################
heatmap.2(as.matrix(pathogens.df[intersect(ipr.genes$symbol, rownames(pathogens.df)),
]),
#betas[head(order(res.lrt$padj), 40),-c(1,2)],
dendrogram = "none",
Colv = F,
colsep = match(c("by.virus.Orsay..Chen.",
                 "Pseudomonas.sp..vs..OP50",
                 "Bacillus.megaterium",
                 "by.M..humicola..N2...24h."),
               N2.experiments.estreses.df$experiment),
sepcolor = "black",
sepwidth = c(0.2, 0.2),
lwd = 1.5,
col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(3),
#col = "redgreen",
trace = "none",
cexRow = 0.6,
srtRow = 20,
srtCol = 40,
#adjCol = 0.5,
cexCol = 0.5,
keysize = 0.4,
ColSideColors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#79AF97FF", "#6A6599FF")[as.numeric(N2.experiments.estreses.df$Category)],
#key.title = "log2FC"
)

venn::venn(list("OrV" = OrV.resp.noSarkies,
                "Bacillus thuriengensis BT247" = Bth.resp.v3,
                "IPR" = ipr.genes$symbol),
           zcolor = pal_jama()(6)[c(1:3,5:6)], #same category colors
           opacity = 0.6,
           #ellipse = T,
           box = F,
           ggplot = T,
           cex = 0.4
) +
  theme(#text = element_text(size = 10),
    axis.title = element_text(size = 22))
Bth.resp.v3
################################################################################
# IMMUNITY
################################################################################
immune.direct.Inferred.degs
heatmap.2(as.matrix(pathogens.df[intersect(immune.direct.Inferred.degs, rownames(pathogens.df)),
]),
#betas[head(order(res.lrt$padj), 40),-c(1,2)],
dendrogram = "none",
Colv = F,
colsep = match(c("by.virus.Orsay..Chen.",
                 "Pseudomonas.sp..vs..OP50",
                 "Bacillus.megaterium",
                 "by.M..humicola..N2...24h."),
               N2.experiments.estreses.df$experiment),
sepcolor = "black",
sepwidth = c(0.2, 0.2),
lwd = 1.5,
col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(3),
#col = "redgreen",
trace = "none",
cexRow = 0.6,
srtRow = 20,
srtCol = 40,
#adjCol = 0.5,
cexCol = 0.5,
keysize = 0.4,
ColSideColors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#79AF97FF", "#6A6599FF")[as.numeric(N2.experiments.estreses.df$Category)],
#key.title = "log2FC"
)
################################################################################
# IMMUNE EFECTORS AFFECTED
################################################################################
# ABC transporters
c("pgp-5", "pmp-2", "pgp-9", "pgp-11", "mrp-5", "pgp-13", "haf-3", "pmp-5", "hmt-1", "haf-9") %in% OrV.resp.noSarkies

################################################################################
# ORV SPECIFIC
################################################################################
orv.specific.v3 <- setdiff(responses.comparison.v3$OrV, unique(unlist(responses.comparison.v3[-1])))

orv.specific.v3 %>%
  writeLines()


heatmap.2(betas.symbol[intersect(orv.specific.v3,
                                 rownames(betas.symbol)),
                       grepl("Groupinfection", colnames(betas.symbol))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "row",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          
          labRow = as.expression(lapply(rownames(betas.symbol[intersect(orv.specific.v3,
                                                                        rownames(betas.symbol)),
                                                              grepl("Groupinfection", colnames(betas.symbol))]),
                                        function(a) bquote(italic(.(a))))),
          cexRow = 0.6,
          adjCol = 0.5,
          srtRow = 1,
          labCol = "",#paste0(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44), " hpi"),
          cexCol = 0.8, 
          srtCol = 20,
          # Add the significance of the log2FC
          cellnote = p.adj.markers.symbol[intersect(orv.specific.v3,
                                                              rownames(betas.symbol)),],
          notecol = "black", 
          notecex = 1,
          #adjCol = 0.5,
          keysize = 0.8,
          key.title = "")

# Distribution of log2FC
hist(betas.symbol[intersect(orv.specific.v3,
                       rownames(betas.symbol)),
             grepl("Groupinfection", colnames(betas.symbol))])

orv.specific.v3.above1 <- apply(abs(betas.symbol[intersect(orv.specific.v3,
                                 rownames(betas.symbol)),
                       grepl("Groupinfection", colnames(betas.symbol))]) > 1,
                       1,
                       any)
orv.specific.v3.above2 <- apply(abs(betas.symbol[intersect(orv.specific.v3,
                                                           rownames(betas.symbol)),
                                                 grepl("Groupinfection", colnames(betas.symbol))]) > 2,
                                1,
                                any)

names(orv.specific.v3.above2[orv.specific.v3.above2 == TRUE]) %>%
  writeLines()

# Only genes with at least a |log2FC| > 2
heatmap.2(betas.symbol[intersect(names(orv.specific.v3.above2[orv.specific.v3.above2 == TRUE]),
                                                        rownames(betas.symbol)),
                       grepl("Groupinfection", colnames(betas.symbol))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "none",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          
          labRow = as.expression(lapply(rownames(betas.symbol[intersect(names(orv.specific.v3.above2[orv.specific.v3.above2 == TRUE]),
                                                                        rownames(betas.symbol)),
                                                              grepl("Groupinfection", colnames(betas.symbol))]),
                                        function(a) bquote(italic(.(a))))),
          cexRow = 0.6,
          adjCol = 0.5,
          srtRow = 1,
          labCol = "",#paste0(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44), " hpi"),
          cexCol = 0.55, 
          srtCol = 20,
          # Add the significance of the log2FC
          cellnote = p.adj.markers.symbol[intersect(names(orv.specific.v3.above2[orv.specific.v3.above2 == TRUE]),
                                                                           rownames(betas.symbol)),
                                                                 ],
          notecol = "black", 
          notecex = 1,
          #adjCol = 0.5,
          keysize = 0.7,
          key.title = "")

heatmap.2(as.matrix(pathogens.df[match(orv.specific.v3,
                                       rownames(pathogens.df)),]),
          #betas[head(order(res.lrt$padj), 40),-c(1,2)],
          dendrogram = "none",
          Colv = F,
          colsep = match(c("by.virus.Orsay..Chen.",
                           "Pseudomonas.sp..vs..OP50",
                           "Bacillus.megaterium",
                           "by.M..humicola..N2...24h."),
                         N2.experiments.estreses.df$experiment),
          sepcolor = "black",
          sepwidth = c(0.2, 0.2),
          #col = c("red", "black", "green"),
          col = c("#c800c9", "white", "#009B9F"),
          #col = "redgreen",
          trace = "none",
          cexRow = 0.7,
          srtRow = 20,
          srtCol = 40,
          #adjCol = 0.5,
          cexCol = 0.5,
          keysize = 0.4,
          ColSideColors = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#79AF97FF", "#6A6599FF"
          )[as.numeric(N2.experiments.estreses.df$Category)],
          #key.title = "log2FC"
)

names(responses.comparison.v3[-1])
