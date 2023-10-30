################################################################################
# Fig 4
# LRT analysis & correlated with virus accumulation
#
# 2023/07/24
# MJ Olmo-Uceda
################################################################################
# A) Venn
################################################################################
# B) Heatmap of the 150 more significant genes
################################################################################
heatmap.2(betas[head(order(res.lrt$padj), 150), grepl("Groupinfection", colnames(betas))],
          #dendrogram = "row",
          dendrogram = "none",
          Colv = F,
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(21),
          #col = "redgreen",
          trace = "none",
          cexRow = 0.5,
          srtRow = 20,
          # with symbol names
          labRow = cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"], 
          labCol = rep("", 11),
          srtCol = 0,
          adjCol = 0.5,
          cexCol = 0,
          
          keysize = 0.5,
          key.title = "log2FC",
          
          # Significativity
          #cellnote = p.adj.markers[head(order(res.lrt$padj), 150),],
          notecol = "black", 
          notecex = 0.8,
)
# C) Clusters present in the 150 genes most significant?
################################################################################
# To which cluster belong?
clusters.allDEGs.lrt.hclust[cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"], ] %>% View()
sort(table(clusters.allDEGs.lrt.hclust[cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"], "Cluster"]),
     decreasing = T)
#  5  2  3  9  4  8  7 11  1 10 14 13 15 16  6 12 17 18 19 
# 52 19 16 16  8  8  5  5  4  3  3  2  2  2  1  1  1  1  1 
# Some of them are in accumulated clusters
sort(table(clusters.allDEGs.lrt.hclust[cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"], "Cluster.manual"]),
     decreasing = T)
# virus-like   late upregulated  from8 upregulated late downregulated            cyclic5 
#         36                 23                 16                 16                  1 
# Reduced clusters df: only the 150 most significant genes
cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"]

clusters.150.mostSig.lrt.hclust <- clusters.allDEGs.lrt.hclust[cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"],]
clusters.150.mostSig.lrt.hclust$cele_id <- cele_trad[match(rownames(clusters.150.mostSig.lrt.hclust), cele_trad$symbol),
                                                     "cele_id"]
## Genes in each interesting cluster
# virus-like
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "virus-like" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames
# [1] "ddn-1"     "eol-1"     "ZC196.3"   "F26F2.1"   "trpl-5"    "pals-5"    "F42C5.3"   "fbxa-182"  "pals-39"  
# [10] "C53A5.9"   "Y39G8B.5"  "droe-8"    "B0507.6"   "F26F2.3"   "pals-4"    "Y46G5A.20" "pals-38"   "pals-37"  
# [19] "T08E11.1"  "F49H6.5"   "lys-3"     "pals-32"   "Y105C5A.9" "F57G4.11"  "pals-26"   "sago-2"    "Y43F8B.12"
# [28] "Y75B8A.39" "C08E3.1"   "pals-27"   "Y6G8.5"    "pals-30"   "pals-14"   "pals-3"    "T07A5.7"   "math-14"  

# from8 upregulated
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "from8 upregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames 
# [1] "math-41"  "glb-1"    "F35E12.6" "F20D6.11" "gst-5"    "ppw-1"    "F55G1.9"  "C35B1.5"  "ZK6.11"   "clec-166"
# [11] "M03A1.3"  "ctsa-4.2" "ugt-19"   "ctsa-4.1" "acly-1"   "math-39" 
# late upregulated
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late upregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames 
# [1] "gst-35"   "asp-14"   "gst-21"   "ugt-21"   "smf-3"    "C32H11.4" "irg-3"    "C18H9.6"  "C32H11.3" "irg-5"   
# [11] "gst-4"    "dod-24"   "cysl-2"   "F55G11.2" "gst-7"    "cyp-29A3" "dct-17"   "clec-143" "nit-1"    "gst-12"  
# [21] "fil-1"    "H20E11.3" "lurp-4"  
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late downregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames()
# [1] "pmp-5"     "lys-4"     "hphd-1"    "T05E12.6"  "asp-13"    "Y39B6A.1"  "lys-6"     "clec-47"   "fat-7"    
# [10] "DH11.2"    "cyp-13A4"  "srh-237"   "Y38H6C.23" "Y38H6C.21" "ftn-1"     "lys-5"    
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster == 5 & !is.na(clusters.150.mostSig.lrt.hclust$Cluster), ] %>% rownames()
# [1] "ahcy-1"    "nhr-68"    "sams-3"    "C23H5.8"   "ugt-62"    "pho-13"    "cpr-4"     "ftn-2"     "cpr-1"    
# [10] "folt-2"    "mthf-1"    "clec-48"   "lbp-6"     "T15B7.1"   "pes-9"     "C09D4.1"   "W01B11.6"  "R09H10.7" 
# [19] "cblc-1"    "ZK1193.2"  "pmt-2"     "kel-8"     "asp-8"     "F54B11.11" "metr-1"    "sptl-1"    "gcs-1"    
# [28] "mtrr-1"    "C30G12.2"  "ent-7"     "F17H10.1"  "nhx-2"     "ech-7"     "ugt-28"    "acs-5"     "C18A11.3" 
# [37] "papl-1"    "pck-2"     "nhr-156"   "fat-6"     "cct-3"     "cdr-4"     "T25B9.9"   "cpr-5"     "T06A1.5"  
# [46] "F58G6.7"   "cth-1"     "mel-32"    "cyp-33E2"  "F58G6.9"   "let-767"   "asp-2"    

# VIRUS-LIKE IN 150 MOST SIGNIFICANT GENES
#"#7D5BA6" (purple)
#"#EFC12E" (mustard) 
#"gray60"
#"#D86357" (dark orange)
#"#AFE198" (lima)
#"#DC9A7E" (pale orange)
infected.color.up <- "#48999D" #"#6EA1CB" #blue
infected.color.down <- "#B52DC2"
ms.virusLike <- ggplot(cluster_lrt[["normalized"]] %>%
         filter(genes %in% clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "virus-like" & 
                                                             !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), "cele_id"]),
       aes(Time,
           value,
           color = Group,
           fill = Group,
           group = Group,
           #shape = Group
       )) +
  # geom_vline(xintercept = as.factor(18),
  #            linetype = 3) +
  # geom_boxplot(aes(group = interaction(Time, Group)),
  #              alpha = 0.1,
  #              ) +
  geom_point(alpha = 0.05) +
  stat_summary(fun = mean,
               geom = "line",
               #aes(linetype = Group)
               ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2.5) +
  # stat_summary(#fun.data = mean_se, 
  #              fun.data = mean_sdl, # standar deviation
  #              fun.args = list(mult = 1), 
  #              geom = "errorbar",
  #              width = 0.2,
  #              alpha = 1) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "ribbon", 
               alpha = 0.15, 
               color = NA, 
               #position = position_dodge(width = 0.5)
               ) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = infected.color.up)) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = infected.color.up)) +
  # scale_shape_manual(values = c("control" = 16,
  #                               "infection" = 1)) +
  # scale_linetype_manual(values = c("control" = "solid",
  #                                  "infection" = "dotted")) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# FROM8 UPREGULATED IN 150 MOST SIGNIFICANT GENES
clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "from8 upregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ]

ms.from8up <- ggplot(cluster_lrt[["normalized"]] %>%
                       filter(genes %in% clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "from8 upregulated" & 
                                                                           !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), "cele_id"]),
                     aes(Time,
                         value,
                         color = Group,
                         fill = Group,
                         group = Group,
                         #shape = Group
                     )) +
  # geom_vline(xintercept = as.factor(18),
  #            linetype = 3) +
  # geom_boxplot(aes(group = interaction(Time, Group)),
  #              alpha = 0.1,
  #              ) +
  geom_point(alpha = 0.05) +
  stat_summary(fun = mean,
               geom = "line",
               #aes(linetype = Group)
               ) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2.5) +
  # stat_summary(fun.data = mean_sdl, # standar deviation
  #              fun.args = list(mult = 1), 
  #              #fun.data = mean_se,  # standar error
  #              geom = "errorbar",
  #              width = 0.2,
  #              alpha = 1) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "ribbon", 
               alpha = 0.15, 
               color = NA, 
               #position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = infected.color.up)) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = infected.color.up)) +
  # scale_shape_manual(values = c("control" = 16,
  #                               "infection" = 1)) +
  # scale_linetype_manual(values = c("control" = "solid",
  #                                  "infection" = "dotted")) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# LATE UPREGULATED IN 150 MOST SIGNIFICANT GENES
ms.lateUp <- ggplot(cluster_lrt[["normalized"]] %>%
         filter(genes %in% clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late upregulated" & 
                                                             !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), "cele_id"]),
       aes(Time,
           value,
           color = Group,
           fill = Group,
           group = Group
       )) +
  # geom_vline(xintercept = as.factor(18),
  #            linetype = 3) +
  # geom_boxplot(aes(group = interaction(Time, Group)),
  #              alpha = 0.1,
  #              ) +
  geom_point(alpha = 0.05) +
  stat_summary(fun = mean,
               geom = "line") +
  stat_summary(fun = mean,
               geom = "point",
               size = 2.5) +
  # stat_summary(#fun.data = mean_se,  
  #              fun.data = mean_sdl, # standar deviation
  #              fun.args = list(mult = 1), 
  #              geom = "errorbar",
  #              width = 0.2,
  #              alpha = 1) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "ribbon", 
               alpha = 0.15, 
               color = NA, 
               #position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = infected.color.up)) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = infected.color.up)) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
# CLUSTER 5 IN 150 MOST SIGNIFICANT GENES
ms.5 <- ggplot(cluster_lrt[["normalized"]] %>%
                 filter(genes %in% clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster == 5, "cele_id"]),
               aes(Time,
                   value,
                   color = Group,
                   fill = Group,
                   group = Group
               )) +
  # geom_vline(xintercept = as.factor(18),
  #            linetype = 3) +
  # geom_boxplot(aes(group = interaction(Time, Group)),
  #              alpha = 0.1,
  #              ) +
  geom_point(alpha = 0.05) +
  stat_summary(fun = mean,
               geom = "line") +
  stat_summary(fun = mean,
               geom = "point",
               size = 2.5) +
  # stat_summary(#fun.data = mean_se,  
  #   fun.data = mean_sdl, # standar deviation
  #   fun.args = list(mult = 1), 
  #   geom = "errorbar",
  #   width = 0.2,
  #   alpha = 1) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "ribbon", 
               alpha = 0.15, 
               color = NA, 
               #position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = infected.color.down)) + #78BE74
  scale_fill_manual(values = c("control" = "black",
                               "infection" = infected.color.down)) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  )
# LATE DOWNREGULATED IN 150 MOST SIGNIFICANT GENES
ms.lateDown <- ggplot(cluster_lrt[["normalized"]] %>%
         filter(genes %in% clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late downregulated" & 
                                                             !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), "cele_id"]),
       aes(Time,
           value,
           color = Group,
           fill = Group,
           group = Group
       )) +
  # geom_vline(xintercept = as.factor(18),
  #            linetype = 3) +
  # geom_boxplot(aes(group = interaction(Time, Group)),
  #              alpha = 0.1,
  #              ) +
  geom_point(alpha = 0.05) +
  stat_summary(fun = mean,
               geom = "line") +
  stat_summary(fun = mean,
               geom = "point",
               size = 2.5) +
  # stat_summary(#fun.data = mean_se,  # standar errors
  #              fun.data = mean_sdl, # standar deviation
  #              fun.args = list(mult = 1), 
  #              geom = "errorbar",
  #              width = 0.2,
  #              alpha = 1) +
  stat_summary(fun.data = mean_sdl, # standar deviation
               fun.args = list(mult = 1), 
               #fun.data = mean_se,  # standar errors
               geom = "ribbon", 
               alpha = 0.15, 
               color = NA, 
               #position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = infected.color.down)) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = infected.color.down)) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank()
        )


wrap_plots(A = ms.virusLike,
           B = ms.from8up,
           C = ms.lateUp,
           D = ms.5,
           E = ms.lateDown,
           design = '
           A
           B
           C
           D
           E'
             )
# Saving the genes forming profiles
lrt.150.clustered.degs <- list("virus-like" = clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "virus-like" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames(),
                               "from8 upregulated" = clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "from8 upregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames(), 
                               "late upregulated" = clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late upregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames(),
                               "late downregulated" = clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster.manual == "late downregulated" & !is.na(clusters.150.mostSig.lrt.hclust$Cluster.manual), ] %>% rownames(),
                               "downregulated" = clusters.150.mostSig.lrt.hclust[clusters.150.mostSig.lrt.hclust$Cluster == 5 & !is.na(clusters.150.mostSig.lrt.hclust$Cluster), ] %>% rownames())
# CORRELATED & ANTICORRELATED GENES
################################################################################
heatmap.2(betas.symbol[intersect(genes_corr_virus, 
                                 rownames(betas.symbol)),
                       grepl("Groupinfection", colnames(betas.symbol))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "none",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          
          labRow = as.expression(lapply(rownames(betas.symbol[intersect(genes_corr_virus, 
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
          cellnote = p.adj.markers.symbol[intersect(genes_corr_virus,
                                                    rownames(betas.symbol)),
          ],
          notecol = "black", 
          notecex = 1,
          #adjCol = 0.5,
          keysize = 0.8,
          key.title = "")

heatmap.2(betas.symbol[intersect(genes_anticorr_virus, 
                                 rownames(betas.symbol)),
                       grepl("Groupinfection", colnames(betas.symbol))], 
          Rowv = T, 
          Colv = FALSE,
          dendrogram = "none",
          col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(100),        
          key = T, 
          trace = "none",
          # symbol but in the same order
          
          labRow = as.expression(lapply(rownames(betas.symbol[intersect(genes_anticorr_virus, 
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
          cellnote = p.adj.markers.symbol[intersect(genes_anticorr_virus,
                                                    rownames(betas.symbol)),
          ],
          notecol = "black", 
          notecex = 1,
          #adjCol = 0.5,
          keysize = 0.8,
          key.title = "")


table(clusters.allDEGs.lrt.hclust[genes_anticorr_virus, "Cluster"])

genes_corr_virus
## Intersections
intersect(genes_corr_virus, ipr.genes$symbol)
# "Y39G8B.5" "ZC196.3" 
intersect(genes_corr_virus, mishra19$symbol) # 0
intersect(genes_corr_virus, zip1.dependent.30min$V1)
# "Y39G8B.5" "T07A5.7" 
intersect(genes_corr_virus, zip1.dependent.4h$V1)
# "C17E4.6" "pes-4"   "taf-6.2" "pxf-1"  
writeLines(intersect(genes_corr_virus, orv.specific.v3), sep = ", ")
# Y65B4A.4, akir-1, unc-57, ztf-6, F46A8.13, W04A8.1, smu-2, ash-2, taf-6.2, pxf-1, ZC196.3
writeLines(intersect(genes_anticorr_virus, orv.specific.v3), sep = ", ")
# tpi-1, F46B6.6
writeLines(intersect(genes_anticorr_virus, common.v3), sep = ", ")



# D) Functional Enrichment of all significant LRT genes
################################################################################
allLRT.wormenrichr.df %>%
  filter(Adjusted.P.value <= 0.05,
         db %in% c("GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GO_Biological_Process_2018",
                   "KEGG_2019")
  ) %>%
  arrange(db, Adjusted.P.value) %>%
  ggplot(.,
         aes(x = fct(ifelse(nchar(gsub(" \\(GO:[0-9]+\\)", "", Term)) > 50,
                            paste(substring(gsub(" \\(GO:[0-9]+\\)", "", Term), 1, 50), "..."),
                            gsub(" \\(GO:[0-9]+\\)", "", Term))),
             # y =  reorder(fct_rev(ifelse(nchar(gsub(" \\(GO:[0-9]+\\)", "", Term)) > 50,
             #                     paste(substring(gsub(" \\(GO:[0-9]+\\)", "", Term), 1, 50), "..."), 
             #                     gsub(" \\(GO:[0-9]+\\)", "", Term))),
             #              Adjusted.P.value),
             y = Adjusted.P.value,
             size = GeneRatio,
             color = db,
             fill = db,
         )) +
  geom_point(shape = 21,
             color = "black") +
  theme_classic() +
  scale_fill_manual(values = db.colors) +
  coord_cartesian(ylim = c(-0.001,0.05)) +
  #scale_x_discrete(expand = expansion(add = c(0.01, 0.01))) +
  #scale_fill_viridis_d(option = "inferno") +
  #limits = c(1,max.Counts)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8,
                                   angle = 10,
                                   color = "black"),
        axis.text.x = element_text(angle = 80,
                                   color = "black",
                                   hjust = 1),
        panel.grid.major.x = element_line(colour = "black",
                                          linewidth = 0.04),
        #panel.spacing.y = unit(0.5, "mm"),
        strip.text = element_blank()
  )

################################################################################
# Functional Enrichment of all LRT genes with enrichGO
################################################################################
# default with clusterProfiler
dotplot(clusterProfiler::simplify(clusterProfiler::enrichGO(as.character(cele_trad[match(all.degs$Profile, cele_trad$cele_id),
                                                                      "geneid"]),
                                                            OrgDb = "org.Ce.eg.db", 
                                                            ont="all", 
                                                            readable=TRUE)),
        #color = geneList.lrt,
        split = "ONTOLOGY") +
  facet_grid(.~ONTOLOGY, 
             scale = "free") +
  theme_classic() +
  coord_flip() +
  scale_color_viridis_c(option = "inferno") +
  theme(#axis.text.y = element_text(size = 8),
    axis.text = element_text(size = 10,
                             angle = 90,
                             hjust = 1),
    legend.position = "top"
  )
## with enrichr (more complete)
## Fixing the legends
### P.adj (fill)
max.padj.filtered <- max(allLRT.wormenrichr.df %>%
                    filter(Adjusted.P.value <= 0.05) %>%
                    dplyr::select(Adjusted.P.value))
min.padj.filtered <- min(allLRT.wormenrichr.df %>%
                           filter(Adjusted.P.value <= 0.05) %>%
                           dplyr::select(Adjusted.P.value))
### Number of genes (size)
max.genes <- max(allLRT.wormenrichr.df %>%
                   filter(Adjusted.P.value <= 0.05) %>%
                   dplyr::select(nGenes))
min.genes <- min(allLRT.wormenrichr.df %>%
                   filter(Adjusted.P.value <= 0.05) %>%
                   dplyr::select(nGenes))
# GO.BP
lrt.bp <- allLRT.wormenrichr.df %>%
  filter(Adjusted.P.value <= 0.05,
         db == "GO_Biological_Process_2018") %>%
  ggplot(.,
         aes(x =  reorder(fct_rev(ifelse(nchar(gsub(" \\(GO:[0-9]+\\)", "", Term)) > 50,
                                 paste(substring(gsub(" \\(GO:[0-9]+\\)", "", Term), 1, 50), "..."), 
                                 gsub(" \\(GO:[0-9]+\\)", "", Term))),
                          GeneRatio),
             y = GeneRatio,
             #size = nGenes,
             fill = Adjusted.P.value,
         )) +
  # geom_tile(color = "black",
  #           height = 0.7
  # ) +
  geom_point(shape = 21,
             size = 4,
             color = "black") +
  theme_classic() +
  #scale_x_discrete(drop = F) +
  # scale_color_manual(values = c("UP" = "#009B9F",
  #                              "DOWN" = "#c800c9")) +
  # scale_fill_manual(values = c("UP" = "#009B9F",
  #                              "DOWN" = "#c800c9")) +
  # scale_fill_gradient(low = "yellow",
  #                     high = "black",
  #                     limits = c(min.padj.filtered,max.padj.filtered)) +
  scale_fill_gradientn(colours = c(viridis::inferno(21)[1],
                                   viridis::inferno(21)[11],
                                   viridis::inferno(21)[21]
                                   ),
                      limits = c(min.padj.filtered,max.padj.filtered)) +
  #limits = c(min.padj.filtered,)) +
  ylim(c(0,1)) +
  #scale_size(range = c(1, 10)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8,
                                   angle = 0),
        axis.text.x = element_text(angle = 80,
                                   hjust = 1),        
        panel.grid.major = element_line(colour = "black",
                                        linewidth = 0.1),
        panel.spacing.y = unit(0.5, "mm")
  )

lrt.mf <- allLRT.wormenrichr.df %>%
  filter(Adjusted.P.value <= 0.05,
         db == "GO_Molecular_Function_2018") %>%
  ggplot(.,
         aes(x =  reorder(fct_rev(ifelse(nchar(gsub(" \\(GO:[0-9]+\\)", "", Term)) > 50,
                                 paste(substring(gsub(" \\(GO:[0-9]+\\)", "", Term), 1, 50), "..."), 
                                 gsub(" \\(GO:[0-9]+\\)", "", Term))),
                          GeneRatio),
             y = GeneRatio,
             #size = GeneRatio,
             fill = Adjusted.P.value,
             #fill = change
         )) +
  # geom_tile(color = "black",
  #           height = 0.7
  # ) +
  geom_point(shape = 21,
             size = 4,
             color = "black") +
  theme_classic() +
  ylim(c(0,1)) +
  #scale_x_discrete(drop = F) +
  scale_fill_gradientn(colours = c(viridis::inferno(21)[1],
                                   viridis::inferno(21)[11],
                                   viridis::inferno(21)[21]
  ),
  limits = c(min.padj.filtered,max.padj.filtered)) +
  # scale_fill_gradient(low = "yellow",
  #                     high = "blue",
  #                     limits = c(min.padj.filtered,max.padj.filtered)) +
  #limits = c(min.padj.filtered,)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 8,
        #                            angle = 10),
        axis.text.x = element_text(angle = 80,
                                   hjust = 1),
        panel.grid.major = element_line(colour = "black",
                                        linewidth = 0.1),
        panel.spacing.y = unit(0.5, "mm")
  )

lrt.cc <- allLRT.wormenrichr.df %>%
  filter(Adjusted.P.value <= 0.05,
         db == "GO_Cellular_Component_2018") %>%
  ggplot(.,
         aes(x =  reorder(fct_rev(ifelse(nchar(gsub(" \\(GO:[0-9]+\\)", "", Term)) > 50,
                                         paste(substring(gsub(" \\(GO:[0-9]+\\)", "", Term), 1, 50), "..."), 
                                         gsub(" \\(GO:[0-9]+\\)", "", Term))),
                          GeneRatio),
             y = GeneRatio,
             #size = GeneRatio,
             fill = Adjusted.P.value,
             #fill = change
         )) +
  # geom_tile(color = "black",
  #           height = 0.7
  # ) +
  geom_point(shape = 21,
             size = 4,
             color = "black") +
  theme_classic() +
  ylim(c(0,1)) +
  #scale_x_discrete(drop = F) +
  scale_fill_gradientn(colours = c(viridis::inferno(21)[1],
                                   viridis::inferno(21)[11],
                                   viridis::inferno(21)[21]
  ),
  limits = c(min.padj.filtered,max.padj.filtered)) +
  # scale_fill_gradient(low = "yellow",
  #                     high = "blue",
  #                     limits = c(min.padj.filtered,max.padj.filtered)) +
  #limits = c(min.padj.filtered,)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 8,
        #                            angle = 10),
        axis.text.x = element_text(angle = 80,
                                   hjust = 1),
        panel.grid.major = element_line(colour = "black",
                                        linewidth = 0.1),
        panel.spacing.y = unit(0.5, "mm")
  )
lrt.kegg <- allLRT.wormenrichr.df %>%
  filter(Adjusted.P.value <= 0.05,
         db == "KEGG_2019") %>%
  ggplot(.,
         aes(x =  reorder(fct_rev(ifelse(nchar(Term) > 50,
                                 paste(substring(Term, 1, 50), "..."), 
                                 Term)),
                          GeneRatio),
             y = GeneRatio,
             #size = GeneRatio,
             fill = Adjusted.P.value,
             #fill = change
         )) +
  geom_point(shape = 21,
             size = 4,
             color = "black") +
  theme_classic() +
  ylim(c(0,1)) +
  scale_fill_gradientn(colours = c(viridis::inferno(21)[1],
                                   viridis::inferno(21)[11],
                                   viridis::inferno(21)[21]
  ),
  limits = c(min.padj.filtered,max.padj.filtered)) +
  # scale_fill_gradient(low = "yellow",
  #                     high = "blue",
  #                     limits = c(min.padj.filtered,max.padj.filtered)) +
  #limits = c(min.padj.filtered,)) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 8,
        #                            angle = 0),
        axis.text.x = element_text(angle = 80,
                                   hjust = 1),
        panel.grid.major = element_line(colour = "black",
                                        linewidth = 0.1),
        panel.spacing.y = unit(0.5, "mm")
  )

(lrt.bp | lrt.mf | lrt.cc | lrt.kegg) 

# horizontal panel
wrap_plots(A = lrt.bp,
           B = lrt.mf,
           C = lrt.cc,
           D = lrt.kegg,
           design = 'AAABBBCDDDD')

wrap_plots(A = lrt.bp,
           B = lrt.mf,
           C = lrt.cc,
           D = lrt.kegg,
           design = 
           'AAABBB
           CCDDDD')

# allLRT.wormenrichr.df %>%
#   filter(Adjusted.P.value <= 0.05,
#          db %in% c(
#          "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GO_Biological_Process_2018")) %>%
#   ggplot(.,
#          aes(x = Term,
#              y = GeneRatio,
#              color = Adjusted.P.value,
#              size = GeneRatio)) +
#   geom_point() +
#   # geom_tile(color = "black",
#   #           height = 1
#   # ) +
#   facet_wrap(~db,
#              scales = "free",
#              ncol = 3) +
#   scale_x_discrete(drop = F) +
#   theme_classic() +
#   coord_flip() +
#   theme(axis.title = element_blank(),
#         #legend.position = "none",
#         axis.text.y = element_text(size = 7,
#                                    angle = 0),
#         axis.text.x = element_blank(),
#         panel.grid.major = element_line(colour = "black",
#                                         linewidth = 0.05),
#         panel.spacing.y = unit(0.6, "mm"))




