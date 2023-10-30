################################################################################
# RESULTS GENE COUNTS WITH STAR
# MJ Olmo-Uceda
# 2023/04/19
################################################################################
# The STAR -geneCounts output is a three column .tab where:
# column 1: gene id
# column 2: unstranded counts
# column 3: stranded-forward counts
# column 4: stranded-reverse counts <-- our libraries

# Loading the data
################################################################################
path.counts <- "data/star.counts/"
count.files <- list.files(path.counts) 

counts <- read.csv(paste0(path.counts, count.files[1]), 
                   header = F, 
                   skip = 4, 
                   sep = "\t")
sample <- unlist(strsplit(count.files[1], split = "_N2_ReadsPerGene.out.tab"))
sample <- paste0(unlist(strsplit(sample, "_")), collapse = ".")

rownames(counts) <- counts$V1
counts.unstranded <- counts
counts.revstranded <- counts

counts.unstranded$V1 <- NULL
counts.unstranded$V3 <- NULL
counts.unstranded$V4 <- NULL
#colnames(counts.unstranded) <- sample

counts.revstranded$V1 <- NULL
counts.revstranded$V2 <- NULL
counts.revstranded$V3 <- NULL
#colnames(counts.revstranded) <- sample

samples <- c(sample)
for (i in 2:length(count.files)){
  sample <- unlist(strsplit(count.files[i], split = "_N2_ReadsPerGene.out.tab"))
  sample <- paste0(unlist(strsplit(sample, "_")), collapse = ".")
  samples[i] <- sample
  c <- read.csv(paste0(path.counts, count.files[i]), 
                header = F, 
                skip = 4, 
                sep = "\t")
  counts.unstranded <- cbind(counts.unstranded, c$V2)
  counts.revstranded <- cbind(counts.revstranded, c$V4)
}
colnames(counts.unstranded) <- samples
colnames(counts.revstranded) <- samples

evolSamples <- c("JUv1580.RL1", "JUv1580.RL2", "JUv2572.10.5", "JUv2572.RL3")
samples %in% evolSamples

# removing evo samples
counts.unstranded <- counts.unstranded %>%
  select(-starts_with("JUv"))
counts.revstranded <- counts.revstranded %>%
  select(-starts_with("JUv"))

# check order
all.equal(colnames(counts.unstranded), infoSamples.complete$Sample)
all.equal(colnames(counts.revstranded), infoSamples.complete$Sample)

dim(counts.unstranded) # 45706    72
dim(counts.revstranded) # 45706    72

# checking the last row is 'transcript_id'
tail(counts.revstranded)
tail(head(counts.revstranded, -1))

# removing the last row 'transcript_ir'
counts.revstranded <- head(counts.revstranded, -1)
counts.unstranded <- head(counts.unstranded, -1)

# changing sample names to "new_names" : {c/i}.{tp}.{replicate}
################################################################################
# We already checked that are in the same order
rownames(infoSamples.complete) <- infoSamples.complete$new_names
colnames(counts.unstranded) <- rownames(infoSamples.complete)
colnames(counts.revstranded) <- rownames(infoSamples.complete)

# write matrix counts of OrVProgression
write.csv(counts.unstranded, 
          file="data/rawCounts_wStarUnstranded_orvProg230419.csv",
          quote = F,
          row.names = T)
write.csv(counts.revstranded, 
          file="data/rawCounts_wStarREVstranded_orvProg230419.csv",
          quote = F,
          row.names = T)

# Sorting the name of the samples per condition and per time
dput(rownames(infoSamples.complete))
sorted.samples <- c(# control
                  "c.0.1", "c.0.2", "c.0.3", 
                  "c.2.1", "c.2.2", "c.2.3", 
                  "c.4.1", "c.4.2", "c.4.3", 
                  "c.6.1", "c.6.2", "c.6.3", 
                  "c.8.1", "c.8.2", "c.8.3", 
                  "c.12.1", "c.12.2", "c.12.3", 
                  "c.18.1", "c.18.2", "c.18.3", 
                  "c.22.1", "c.22.2", "c.22.3", 
                  "c.26.1", "c.26.2", "c.26.3", 
                  "c.32.1", "c.32.2", "c.32.3", 
                  "c.38.1", "c.38.2", "c.38.3", 
                  "c.44.1", "c.44.2", "c.44.3", 
                  # infected
                  "i.0.1", "i.0.2", "i.0.3", 
                  "i.2.1", "i.2.2", "i.2.3", 
                  "i.4.1", "i.4.2", "i.4.3", 
                  "i.6.1", "i.6.2", "i.6.3", 
                  "i.8.1", "i.8.2", "i.8.3",
                  "i.12.1", "i.12.2", "i.12.3", 
                  "i.18.1", "i.18.2", "i.18.3", 
                  "i.22.1", "i.22.2", "i.22.3", 
                  "i.26.1", "i.26.2", "i.26.3",
                  "i.32.1", "i.32.2", "i.32.3", 
                  "i.38.1", "i.38.2", "i.38.3", 
                  "i.44.1", "i.44.2", "i.44.3"
                  )

counts.revstranded.sorted <- counts.revstranded[,sorted.samples]
write.csv(counts.revstranded.sorted, 
          file="data/rawCounts.sorted_wStarREVstranded_orvProg230420.csv",
          quote = F,
          row.names = T)
# We will use for th next analysis only the revstranded data

# Low count filter:
## Filtering genes with expressions of min 1 CPM in at least 3 samples
################################################################################
dim(counts.unstranded[rowSums(counts.unstranded == 0) == ncol(counts.unstranded),]) # 13611
dim(counts.revstranded[rowSums(counts.revstranded == 0) == ncol(counts.revstranded),]) # 11070
dim(counts.revstranded.sorted[rowSums(counts.revstranded.sorted == 0) == ncol(counts.revstranded.sorted),]) # 11070

cts.filtered <- counts.revstranded.sorted[rowSums(edgeR::cpm(counts.revstranded.sorted) >= 1) > 2,]
dim(cts.filtered) # 16933    72

write.csv(cts.filtered, 
          file="data/filteredCounts1CPM_wStarREVstranded_orvProg230420.csv",
          quote = F,
          row.names = T)

dput(rownames(infoSamples.complete))
# sort infoSamples.complete by condition and time
infoSamples.complete <- infoSamples.complete[sorted.samples,]

# checking order
all.equal(colnames(cts.filtered), rownames(infoSamples.complete))

colSums(counts.revstranded)
