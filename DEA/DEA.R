################################################################################
## Time Course experiment
# Following the suggested in the DESeq2 vignette:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
################################################################################
dds.lrt <- DESeqDataSetFromMatrix(cts.filtered,
                                  infoSamples.complete,
                                  design = ~ Group + Time + Group:Time)
dds.lrt <- DESeq(dds.lrt,
                 reduced = ~ Group + Time,
                 test = "LRT")

resultsNames(dds.lrt)
plotDispEsts(dds.lrt)

res.lrt <- results(dds.lrt)
# sig results
sig.res.lrt <- res.lrt %>%
  as.data.frame() %>%
  filter(padj <= 0.05) %>%
  arrange(padj) 
sig.res.lrt <- merge(sig.res.lrt,
                     cele_trad,
                     by.x = 0,
                     by.y = "cele_id",
                     sort = F,
                     all.x = T)
# write.results.profile
write.xlsx(sig.res.lrt,
           "results/DEgenes/sig.lrt.xlsx")
# results filtering genes by padj
### padj <= 0.01
results(dds.lrt)  %>%
  as.data.frame() %>%
  filter(padj <= 0.01) %>%
  arrange(padj) %>%
  rownames() %>%
  write.table(.,
              file = "results/DEgenes/DEgenes.lrt.padj0.01.txt",
              quote = F,
              row.names = F,
              col.names = F)

### padj <= 0.05, sorted by significativity
genes.lrt <- results(dds.lrt)  %>%
  as.data.frame() %>%
  filter(padj <= 0.05) %>%
  arrange(padj) %>%
  rownames()

results(dds.lrt)  %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  rownames() %>%
  write.table(.,
              file = "results/DEgenes/DEgenes.lrt.padj0.05.txt",
              quote = F,
              row.names = F,
              col.names = F)

pdf("results/DEgenes/deg.profiles/deg.genes.lrt.padj0.05.pdf"
    #width = 8.3,
    #height = 5.8
    )
# Plot the profiles of the significant DEGs
for (gene in genes.lrt){
  print(plotCounts(dds.lrt,
             gene, 
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
    ggtitle(gene) +
    ylab("Expression level (log10(DESeq2 normalized counts))") +
    xlab("Time (hpi)") +
    scale_color_manual(values = c("control" = "black",
                                  "infection" = "red")) +
    scale_fill_manual(values = c("control" = "black",
                                 "infection" = "red")) +
    scale_y_log10() +
    theme_classic())
}
dev.off()

################################################################################
# Wald tests for the log2 fold changes at individual time points
################################################################################
resultsNames(dds.lrt)
# [1] "Intercept"                  "Group_infection_vs_control" "Time_2_vs_0"               
# [4] "Time_4_vs_0"                "Time_6_vs_0"                "Time_8_vs_0"               
# [7] "Time_12_vs_0"               "Time_18_vs_0"               "Time_22_vs_0"              
# [10] "Time_26_vs_0"               "Time_32_vs_0"               "Time_38_vs_0"              
# [13] "Time_44_vs_0"               "Groupinfection.Time2"       "Groupinfection.Time4"      
# [16] "Groupinfection.Time6"       "Groupinfection.Time8"       "Groupinfection.Time12"     
# [19] "Groupinfection.Time18"      "Groupinfection.Time22"      "Groupinfection.Time26"     
# [22] "Groupinfection.Time32"      "Groupinfection.Time38"      "Groupinfection.Time44"     

# At x hpi
results(dds.lrt,
        name = "Groupinfection.Time18",
        test = "Wald") %>%
  as.data.frame() %>%
  filter(padj <= 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  dplyr::select(c(log2FoldChange)) %>%
  write.table(quote = F,
              row.names = T)

# plot interesting gene
g = "CELE_F58B3.1"
plotCounts(dds.lrt,
           g, 
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
  ggtitle(g) +
  ylab("Expression level (log10(DESeq2 normalized counts))") +
  xlab("Time (hpi)") +
  scale_color_manual(values = c("control" = "black",
                                "infection" = "red")) +
  scale_fill_manual(values = c("control" = "black",
                               "infection" = "red")) +
  scale_y_log10() +
  theme_classic()

################################################################################
## Saving all results by timepoint in the list --> res.tp
res.tp <- c()
for (res in dput(resultsNames(dds.lrt)[grepl("Groupinfection", resultsNames(dds.lrt))])){
  print(res)
  res. <- results(dds.lrt,
                            name = paste0(res),
                            test = "Wald") %>%
    as.data.frame() 
  res.tp[[res]] <- merge(res., cele_trad,
        by.x = 0,
        by.y = "cele_id",
        sort = F,
        all.x = T)
  sig.results.degs.tp[[res]] <- res.
}

dput(names(res.tp))
# c("Groupinfection.Time2", "Groupinfection.Time4", "Groupinfection.Time6", 
#   "Groupinfection.Time8", "Groupinfection.Time12", "Groupinfection.Time18", 
#   "Groupinfection.Time22", "Groupinfection.Time26", "Groupinfection.Time32", 
#   "Groupinfection.Time38", "Groupinfection.Time44")

names(res.tp) <- c("2 hpi", "4 hpi", "6 hpi", 
                   "8 hpi", "12 hpi", "18 hpi", 
                   "22 hpi", "26 hpi", "32 hpi", 
                   "38 hpi", "44 hpi")
##############################################################################
## Saving all the significant (padj <= 0.05) DEGs in separate sheets of a excel
wb <- createWorkbook()
for (res in dput(resultsNames(dds.lrt))){
  addWorksheet(wb, 
               paste0(res))
  writeData(wb, 
            paste0(res),
            results(dds.lrt,
                    name = paste0(res),
                    test = "Wald") %>%
              as.data.frame() %>%
              filter(padj <= 0.05) %>%
              arrange(desc(abs(log2FoldChange))),
            row.names = T)
}
saveWorkbook(wb, "results/DEgenes/degs.lrt.byTime.xlsx", overwrite = TRUE)

## Saving significant genes per time
degs.tp <- c()
for (res in dput(resultsNames(dds.lrt)[grepl("Groupinfection", resultsNames(dds.lrt))])){
  print(res)
  degs.tp[[res]] <- results(dds.lrt,
                            name = paste0(res),
                            test = "Wald") %>%
    as.data.frame() %>%
    filter(padj <= 0.05) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    rownames()
}

# list of significant results by time (Wald test)
sig.results.degs.tp <- c()
# list of sig gene names
degs.tp <- c() 
degs.tp.cele <- c()
for (res in dput(resultsNames(dds.lrt)[grepl("Groupinfection", resultsNames(dds.lrt))])){
  print(res)
  res. <- results(dds.lrt,
                  name = paste0(res),
                  test = "Wald") %>%
  as.data.frame() %>%
  filter(padj <= 0.05) %>%
  arrange(desc(abs(log2FoldChange)))
  res. <- merge(res., cele_trad,
                by.x = 0,
                by.y = "cele_id",
                sort = F,
                all.x = T)
  sig.results.degs.tp[[res]] <- res.
  degs.tp[[res]] <- res.$geneid
  degs.tp.cele[[res]] <- res.$Row.names
}

## Split in UP/DOWN
degs.tp.UP <- c()
degs.tp.DOWN <- c()
for (res in names(sig.results.degs.tp)){
  print(res)
  degs.tp.UP[[res]] <- sig.results.degs.tp[[res]][sig.results.degs.tp[[res]]$log2FoldChange > 0, "geneid"]
  degs.tp.DOWN[[res]] <- sig.results.degs.tp[[res]][sig.results.degs.tp[[res]]$log2FoldChange < 0, "geneid"]
}

degs.tp.symbol <- list()
for (tp in names(degs.tp.cele)){
  degs.tp.symbol[[tp]] <- cele_trad[match(degs.tp.cele[[tp]], cele_trad$cele_id),
                                    "symbol"]
}

wb.s <- createWorkbook()
for (tp in names(degs.tp.symbol)){
  addWorksheet(wb.s, 
               paste0(tp))
  writeData(wb.s, 
            paste0(tp),
            degs.tp.symbol[[tp]],
            row.names = T)
}

saveWorkbook(wb.s, "results/DEgenes/degs.symbol.byTime.xlsx", overwrite = TRUE)

################################################################################
# Matrix with the log2FC at each time (even if they are not significant in a 
# specific timepoint)
################################################################################
betas <- coef(dds.lrt)
heatmap.2(betas[head(order(res.lrt$padj), 100), grepl("Groupinfection", colnames(betas))],
          #betas[head(order(res.lrt$padj), 40),-c(1,2)],
          dendrogram = "row",
          Colv = F,
          col = colorRampPalette(c("magenta", "white", "blue"))(21),
          # col = c("#C75DAA", "#CC6EB1", "#D07DB7", "#D58CBE", "#D99BC5", "#DEA9CC", 
          #         "#E2B8D3", "#E6C6DB", "#E9D4E2", "#EDE2E9", "#FFFFFF", "#DBE8E8", 
          #         "#C6DFE0", "#AFD6D7", "#97CDCF", "#7DC5C7", "#5EBCBF", "#31B4B7", 
          #         "#00ACAF", "#00A3A7", "#009B9F"),
          #col = "redgreen",
          trace = "none",
          cexRow = 0.6,
          labRow = cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 100),]), cele_trad$cele_id), "symbol"],
          labCol = rep(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44),2),
          srtCol = 0,
          adjCol = 0.5,
          keysize = 0.7,
          key.title = "log2FC")

################################################################################
# Matrix with the log2FC 
################################################################################
heatmap.2(betas[head(order(res.lrt$padj), 150), 
                grepl("Groupinfection", colnames(betas))],
  #betas[head(order(res.lrt$padj), 40),-c(1,2)],
  dendrogram = "row",
  Colv = F,
  #col = colorRampPalette(c("#D81B60", "white", "#004D40"))(21),
  #col = colorRampPalette(c("#64C28D", "white", "#AF00DE"))(21),
  #col = colorRampPalette(c("#C75DAA", "white", "#009B9F"))(21),
  #col = colorRampPalette(c("#DB00DD", "white", "#1CDFDC"))(21),
  col = colorRampPalette(c("#c800c9", "white", "#009B9F"))(21),
  #col = "redgreen",
  trace = "none",
  cexRow = 0.4,
  srtRow = 20,
  # with symbol names
  labRow = cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 150),]), cele_trad$cele_id), "symbol"], 
  #labCol = rep(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44),2),
  labCol = rep("", 11),
  srtCol = 0,
  adjCol = 0.5,
  cexCol = 0,
  keysize = 0.7,
  key.title = "log2FC")

heatmap.bp(betas[head(order(res.lrt$padj), 50), grepl("Groupinfection", colnames(betas))],
           rlabels = T)

betas[head(order(res.lrt$padj), 100), grepl("Groupinfection", colnames(betas))] %>%
  rownames() %>%
  writeLines()
# stringdb: https://string-db.org/cgi/network?taskId=bIrEN1hjQ2gM&sessionId=b3ljEKslu3oT

cele_trad[match(rownames(res.lrt[head(order(res.lrt$padj), 100),]), cele_trad$cele_id), "symbol"] %>%
  writeLines()

write.csv(betas[,grepl("Groupinfection", colnames(betas))],
          file = "results/DEgenes/betas.lrt.csv",
          quote = F,
          row.names = T)

################################################################################
# Extracting number of significant genes by timepoint
################################################################################
res.name <- c()
n.sDEGs <- c() 
up <- c()
up.lfc1 <- c()
up.lfc2 <- c()
up.lfc3 <- c()
down <- c()
down.lfc1 <- c()
down.lfc2 <- c()
down.lfc3 <- c()
for (res in resultsNames(dds.lrt)[grepl("Groupinfection", resultsNames(dds.lrt))]){
  print(res)
  res.name <- c(res.name, res)
  res. = results(dds.lrt,
                 name = paste0(res),
                 test = "Wald") %>%
    as.data.frame() %>%
    filter(padj <= 0.05)
  n.sDEGs <- c(n.sDEGs, nrow(res.))
  up <- c(up, nrow(res.[res.$log2FoldChange > 0,]))
  up.lfc1 <- c(up.lfc1, nrow(res.[res.$log2FoldChange > 1,]))
  up.lfc2 <- c(up.lfc2, nrow(res.[res.$log2FoldChange > 2,]))
  up.lfc3 <- c(up.lfc3, nrow(res.[res.$log2FoldChange > 3,]))
  down <- c(down, nrow(res.[res.$log2FoldChange < 0,]))
  down.lfc1 <- c(down.lfc1, nrow(res.[res.$log2FoldChange < -1,]))
  down.lfc2 <- c(down.lfc2, nrow(res.[res.$log2FoldChange < -2,]))
  down.lfc3 <- c(down.lfc3, nrow(res.[res.$log2FoldChange < -3,]))
}

nDEGs <- data.frame(results = res.name,
                    time = factor(c(2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44)),
                    n.sDEGs = n.sDEGs,
                    UP = up,
                    UP.lfc1 = up.lfc1,
                    UP.lfc2 = up.lfc2,
                    UP.lfc3 = up.lfc3,
                    DOWN = down,
                    DOWN.lfc1 = down.lfc1,
                    DOWN.lfc2 = down.lfc2,
                    DOWN.lfc3 = down.lfc3
                    )

# betas with symbol as rownames
betas.symbol <- betas[,colnames(betas)[grepl("Groupinfection.", colnames(betas))]]
rownames(betas.symbol) <- ifelse(!is.na(cele_trad[match(rownames(betas.symbol), cele_trad$cele_id),
                                        "symbol"]),
                              cele_trad[match(rownames(betas.symbol), cele_trad$cele_id),
                                        "symbol"],
                              rownames(betas.symbol))

betas["CELE_F52B11.4", colnames(betas)[grepl("Groupinfection.", colnames(betas))]]


for (tp in names(res.tp)){
  print(tp)
  print(res.tp[[paste0(tp)]] %>%
    filter(padj < 0.05) %>%
    pull(padj) %>%
    max())
}
