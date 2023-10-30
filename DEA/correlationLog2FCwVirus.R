################################################################################
# correlationlog2FCwVirus.R
#
# from log2FC [betas matrix](filtered counts genes with 1CPM in at least 3 samples)
# Testing some correlations between the expression of the worm genes and virus
#
# MJ Olmo-Uceda
# 2023/06/27
################################################################################
## Formating the virus data
str(infoSamples.complete)
# Select the data to use as representative of viral concentration in the sample
# taking into account that we will use the log2FC(geneInf/geneCtr) of the times
# 2 to 44
virus.transcript.summary <- infoSamples.complete %>%
  filter(Group == "infection",
         Time %in% c(2,4,6,8,12,18,22,26,32,38,44)) %>%
  group_by(Time) %>%
  summarize(mean.virus_mapped = mean(mapped_viral_reads_percent),
            sd.virus_mapped = sd(mapped_viral_reads_percent),
            mean.RNA1 = mean(RNA1_OrV_numReads/TotalReads_clean_2p150b),
            sd.RNA1 = sd(RNA1_OrV_numReads/TotalReads_clean_2p150b),
            mean.RNA2 = mean(RNA2_OrV_numReads/TotalReads_clean_2p150b),
            sd.RNA2 = sd(RNA2_OrV_numReads/TotalReads_clean_2p150b),
            mean.ratioRNA2.RNA1.sizeRelat = mean(ratioRNA2.RNA1.sizeRelat),
            sd.ratioRNA2.RNA1.sizeRelat = sd(ratioRNA2.RNA1.sizeRelat))

str(virus.transcript.summary)

# Overall profile of the virus accumulation
ggplot() +
  geom_point(data = virus.transcript.summary,
             aes(x = as.numeric(as.character(Time)),
                 y = mean.virus_mapped)) +
  geom_line(data = virus.transcript.summary,
            aes(x = as.numeric(as.character(Time)),
                y = mean.virus_mapped)) +
  geom_ribbon(data = virus.transcript.summary,
              aes(x = as.numeric(as.character(Time)),
                  ymin = ifelse(mean.virus_mapped - sd.virus_mapped < 0,
                                0, mean.virus_mapped - sd.virus_mapped),
                  ymax = mean.virus_mapped + sd.virus_mapped),
              alpha = 0.2,
              fill = "black") +
  xlab("Time (hpi)") +
  ylab("% virus in sample (RNAseq)") +
  theme_classic()

################################################################################
# The correlation is calculated with all the genes but finally we only will consider
# the significant DEGs 
cor.virus.lfc <- c()
cor.virus.lfc.pvalue <- c()

cor.RNA1.lfc <- c()
cor.RNA1.lfc.pvalue <- c()

cor.RNA2.lfc <- c()
cor.RNA2.lfc.pvalue <- c()

cor.ratio.lfc <- c()
cor.ratio.lfc.pvalue <- c()

for (gene in rownames(betas.symbol)){
  ## with the %mapped to virus
  temp <- cor.test(virus.transcript.summary$mean.virus_mapped,
                                    unlist(betas.symbol[gene,]),
                                    #method = "spearman"
                                    )
  cor.virus.lfc[[gene]] <- temp$estimate
  cor.virus.lfc.pvalue[[gene]] <- temp$p.value
  
  ## specific to vRNA1
  temp <- cor.test(virus.transcript.summary$mean.RNA1,
                   unlist(betas.symbol[gene,]),
                   #method = "spearman"
                   )
  cor.RNA1.lfc[[gene]] <- temp$estimate
  cor.RNA1.lfc.pvalue[[gene]] <- temp$p.value
  
  ## specific to vRNA2
  temp <- cor.test(virus.transcript.summary$mean.RNA2,
                   unlist(betas.symbol[gene,]),
                   #method = "spearman"
  )
  cor.RNA2.lfc[[gene]] <- temp$estimate
  cor.RNA2.lfc.pvalue[[gene]] <- temp$p.value
  
  ## specific to ratioRNA2/RNA1
  temp <- cor.test(virus.transcript.summary$mean.ratioRNA2.RNA1.sizeRelat,
                   unlist(betas.symbol[gene,]),
                   #method = "spearman"
  )
  cor.ratio.lfc[[gene]] <- temp$estimate
  cor.ratio.lfc.pvalue[[gene]] <- temp$p.value
}
# Generate the df
cor.lfc <- data.frame(
  w.virus = unlist(cor.virus.lfc),
  w.virus.pval = unlist(cor.virus.lfc.pvalue),
  w.virus.FDR = p.adjust(unlist(cor.virus.lfc.pvalue),
                         method = "BH"),
  w.RNA1 = unlist(cor.RNA1.lfc),
  w.RNA1.pval = unlist(cor.RNA1.lfc.pvalue),
  w.RNA1.FDR = p.adjust(unlist(cor.RNA1.lfc.pvalue),
                        method = "BH"),
  
  w.RNA2 = unlist(cor.RNA2.lfc),
  w.RNA2.pval = unlist(cor.RNA2.lfc.pvalue),
  w.RNA2.FDR = p.adjust(unlist(cor.RNA2.lfc.pvalue),
                        method = "BH"),
  w.ratioRNA2.RNA1 = unlist(cor.ratio.lfc),
  w.ratioRNA2.RNA1.pval = unlist(cor.ratio.lfc.pvalue),
  w.ratioRNA2.RNA1.FDR =  p.adjust(unlist(cor.ratio.lfc.pvalue),
                                   method = "BH")
)

plot(density(p.adjust(cor.lfc$w.virus.pval)))
plot(density(cor.lfc$w.virus.pval))

# Remove ".cor" from gene name
rownames(cor.lfc) <- str_remove_all(rownames(cor.lfc),
                                    ".cor")
cor.test(virus.transcript.summary$mean.RNA2,
         unlist(betas.symbol[gene,]))

hist(cor.lfc$w.virus.FDR,
    breaks = 500)
hist(apply(cor.lfc[,c("w.RNA1","w.RNA2")],
           1, 
           sd),
     breaks = 350)

min(cor.lfc$w.ratioRNA2.RNA1.FDR)
plot(cor.lfc$w.virus,
     cor.lfc$w.virus.FDR,
     pch = ".")
abline(h = 0.01)

cor.lfc[cor.lfc$w.virus.FDR <= 0.16,]

cor.lfc[abs(cor.lfc$w.virus) >= 0.9 & cor.lfc$w.virus.FDR <= 0.05,] %>%
  nrow()

plot(x = betas.symbol["daf-2",],
     y = virus.transcript.summary$mean.virus_mapped)

virus.transcript.summary[which.max(abs(virus.transcript.summary$mean.RNA1 - virus.transcript.summary$mean.ratioRNA2.RNA1.sizeRelat)),]

heatmap.2(as.matrix(dist(cor.lfc$w.virus)),
          col = "redgreen")

# Positive correlated
rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]) %>%
  #length()
  writeLines()
# Negative correlated
rownames(cor.lfc[cor.lfc$w.virus <= -0.8,]) %>%
  #length()
  writeLines()
betas.symbol["daf-2",]
p.adj.markers.symbol["daf-2",]

################################################################################
# Plot the genes correlated above threshold
################################################################################
## Anticorrelated
ggplot(cluster_lrt[["normalized"]] %>%
             filter(genes %in% cele_trad[match(rownames(cor.lfc[cor.lfc$w.virus <= -0.8,]), 
                                               cele_trad$symbol), 
                                         "cele_id"]
             ),
             aes(Time,
                 value,
                 color = Group,
                 fill = Group,
                 group = Group
             )) +
        geom_jitter(alpha = 0.1) +
        stat_summary(fun = mean,
                     geom = "line",
                     alpha = 0.8) +
        stat_summary(fun = mean,
                     geom = "point",
                     size = 3,
                     alpha = 0.8) +
        stat_summary(fun.data = mean_se,
                     geom = "errorbar",
                     width = 0.2,
                     alpha = 1) +
        scale_color_manual(values = c("control" = "black",
                                      "infection" = "#7D5BA6")) +
        scale_fill_manual(values = c("control" = "black",
                                     "infection" = "#7D5BA6")) +
        ylab("Z-score") +
        xlab("hpi") +
        #ggtitle(paste0("Cluster ", component)) +
        theme_classic() +
        theme(legend.position = "none")
## Correlated
ggplot(cluster_lrt[["normalized"]] %>%
         filter(genes %in% cele_trad[match(rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]), 
                                                                          cele_trad$symbol), 
                                                                    "cele_id"]
         ),
       aes(Time,
           value,
           color = Group,
           fill = Group,
           group = Group
       )) +
  geom_point(#position = "jitter"
    alpha = 0.1) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               alpha = 0.8,
               #shape = "-",
               size = 3) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 0.3,
               alpha = 0.8) +
  scale_color_manual(values = c("control" = "black",
                                "infection" = "#7D5BA6")) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = "#7D5BA6")) +
  ylab("Z-score") +
  xlab("hpi") +
  theme_classic() +
  theme(legend.position = "none")

apply(betas.symbol[rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]),], 2, scale)

# Correlated
rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]) %>%
  #writeLines()
  length() # 223

# Anticorrelated
rownames(cor.lfc[cor.lfc$w.virus <= -0.8,]) %>%
  #writeLines()
  length() # 78

# Correlated
rownames(cor.lfc[cor.lfc$w.virus >= 0.85,]) %>%
  #writeLines()
  length() # 60

# Anticorrelated
rownames(cor.lfc[cor.lfc$w.virus <= -0.85,]) %>%
  #writeLines()
  length() # 27

intersect(cor.df$symbol, rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]))
################################################################################
# Plot but log2FoldChange
################################################################################
betas.zscore.symbol <- scale(betas.symbol)
colnames(betas.zscore.symbol) <- paste("hpi.", c(2,4,6,8,12,18,22,26,32,38,44),
                                       sep = "")
betas.zscore.summary.symbol

betas.zscore.symbol[rownames(cor.lfc[cor.lfc$w.virus >= 0.8,]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.zscore.symbol),
               values_to = "z.score",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "hpi.")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = z.score,
             #group = gene
             )) +
  geom_point(alpha = 0.01) +
  ylim(c(-1.5,1.5)) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.3,
               alpha = 0.8) +
  theme_classic()

betas.symbol[rownames(cor.lfc[cor.lfc$w.virus >= 0.80,]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  theme_classic()

betas.symbol[rownames(cor.lfc[cor.lfc$w.virus <= - 0.8 & 
                                rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                                       cele_trad$cele_id), 
                                                                 "symbol"],]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC,
             #group = gene
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  theme_classic()

betas.symbol[rownames(cor.lfc[cor.lfc$w.virus <= - 0.8,]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC,
             #group = gene
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  theme_classic()

# IQR  
# unlist(apply(betas.zscore.symbol, 2, function(x) quantile(x, c(0.25, 0.75)))[1,])
################################################################################
# Taking into account only the genes that were identified as DEGs (padj <= 0.05)
################################################################################
# all significant DEGs
length(unique(unlist(all.degs))) # 2498
## All the sig DEGs has symbol
table(is.na(cele_trad[match(unique(unlist(all.degs)), cele_trad$cele_id), "symbol"]))

## DEGs Correlated with virus accumulation
betas.symbol[rownames(cor.lfc[cor.lfc$w.virus >=  0.8 & 
                                rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                                       cele_trad$cele_id), 
                                                                 "symbol"],]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC,
             #group = gene
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  theme_classic()

## Anticorrelated
betas.symbol[rownames(cor.lfc[cor.lfc$w.virus <= -0.8 & 
                                rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                                       cele_trad$cele_id), 
                                                                 "symbol"],])
             ,] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC,
             #group = gene
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  # geom_ribbon(aes(x = c(2,4,6,8,12,18,22,26,32,38,44),
  #                 ymin = as.vector(apply(betas.zscore.symbol, 2, function(x) quantile(x, c(0.25, 0.75)))[1,]),
  #                 ymax = as.vector(apply(betas.zscore.symbol, 2, function(x) quantile(x, c(0.25, 0.75)))[2,])),
  #             fill = "black",
  #             alpha = 0.1) +
  theme_classic()

## DEGs Correlated with virus accumulation
betas.symbol[rownames(cor.lfc[cor.lfc$w.virus >=  0.8 & 
                                rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                                       cele_trad$cele_id), 
                                                                 "symbol"],]),] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(.,
               cols = colnames(betas.symbol),
               values_to = "log2FC",
               names_to = "time") %>%
  mutate(Time = str_remove_all(time, "Groupinfection.Time")) %>%
  ggplot(.,
         aes(x = as.numeric(Time),
             y = log2FC,
             #group = gene
         )) +
  geom_line(aes(group = gene),
            alpha = 0.05) +
  geom_point(alpha = 0.01) +
  stat_summary(fun = mean,
               geom = "line",
               alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 2,
               alpha = 0.8) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar",
               width = 1,
               alpha = 0.8) +
  theme_classic()

## Which are?
## Correlated: 45
rownames(cor.lfc[cor.lfc$w.virus >= 0.9 & 
                   rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                          cele_trad$cele_id), 
                                                    "symbol"],]) %>%
  writeLines()
  length()
## Anticorrelated: 13
rownames(cor.lfc[cor.lfc$w.virus <=  -0.8 & 
                   rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                          cele_trad$cele_id), 
                                                    "symbol"],]) %>%
  writeLines()
  length()


## This version1 with only DEGs seems the most logical
rownames(cor.lfc[cor.lfc$w.virus >=  0.8 & 
                   rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                          cele_trad$cele_id), 
                                                    "symbol"],])

"clic-1" %in% rownames(cor.lfc[cor.lfc$w.virus >=  0.8 & 
                                 rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                                        cele_trad$cele_id), 
                                                                  "symbol"],])
rownames(cor.lfc[cor.lfc$w.virus >=  0.8,]) %>%
  writeLines()
# correlation & DEGs
rownames(cor.lfc[cor.lfc$w.virus >=  0.8 & 
                   rownames(cor.lfc) %in% cele_trad[match(unique(unlist(all.degs)), 
                                                          cele_trad$cele_id), 
                                                    "symbol"],]) %>%
  writeLines()

# 
# clathrin-mediated endocytosis (act-5, apb-1, chc-1, dyn-1, hsp-1, rab-11.1 and wsp-1) 
"tol-1" %in% rownames(cor.lfc[abs(cor.lfc$w.virus) >= 0.8,])
"dcr-1" %in% rownames(cor.lfc[cor.lfc$w.ratioRNA2.RNA1 >= 0.5,])

# is a gene significantly DE
"tol-1" %in% cele_trad[match(unique(unlist(all.degs)), cele_trad$cele_id), "symbol"]

cor.lfc["tol-1",]

plotCounts(dds.lrt,
           cele_trad[cele_trad$symbol == "pals-5", "cele_id"], 
           intgroup = c("Time","Group"),
           returnData = TRUE) %>%
  ggplot(., 
         aes(x = as.numeric(as.character(Time)),
             y = count, 
             color = Group, 
             group = Group)) + 
  geom_point(alpha = 0.7) + 
  stat_summary(fun = mean, 
               geom = "line") +
  #ggtitle(gene) +
  ylab("Expression level (log10(DESeq2 normalized counts))") +
  xlab("Time (hpi)") +
  scale_color_manual(values = c("control" = "black",
                                "infection" = "red")) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = "red")) +
  scale_y_log10() +
  theme_classic()

cor.df[cor.df$virus >= 0.8 & !is.na(cor.df$virus), "symbol"] %>%
  writeLines()

