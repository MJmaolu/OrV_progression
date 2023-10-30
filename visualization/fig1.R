# fig1.R
# MJ Olmo-Uceda
# 2023/04/18
################################################################################
# B) qPCR rna2 + Viral fragments (%) per sample
ggplot(vRNA.load,
       aes(x = Time,
           y = Viral.load)) +
  geom_point(#size = 1
             ) +
  geom_point(aes(y = mapped_viral_reads_percent * 5 * 10 ** 6),
             shape = 1,
             #fill = "#1D7874",
             #color = "#5F0F40", 
             color = "gray20",
             #color = "yellow",
             #size = 1
             ) +
  # mean values (3 rep) as lines
  geom_line(data = summary.vRNA.load,
            aes(x = Time,
                y = mean.viralLoad)) +
  geom_line(data = summary.vRNA.load[!is.na(summary.vRNA.load$mean.fragments),],
            aes(x = Time,
                y = (mean.fragments) * 5 * 10 ** 6),
            #color = "#1D7874",
            color = "gray20",
            #color = "yellow",
            #linewidth = 1.4,
            linetype = 2) +
  xlab("Time (hpi)") +
  ylab("Viral load") +
  theme_classic() +
  theme(axis.text =  element_blank(),
        axis.title = element_blank()
        )
####Con y-axis con valores transcriptómica a la derecha
ggplot(vRNA.load %>% filter(Time != 1), 
       aes(x = Time)) +
  geom_point(aes(y = Viral.load), 
             color = "black") +
  # Add a secondary y-axis
  scale_y_continuous(
    name = "",
    breaks = scales::breaks_pretty(),
    labels = scales::label_number(scale = 10^-7),
    sec.axis = sec_axis(
      ~ . / (5 * 10^6),  # Scale Variable2 by dividing by 5 * 10^6
      name = "",
      breaks = scales::breaks_pretty(),
      labels = scales::label_number()
    )) +
  # Add white points for Variable2
  geom_point(aes(y = mapped_viral_reads_percent * 5 * 10 ** 6),
             shape = 1) +
  # mean values (3 rep) as lines
  geom_line(data = summary.vRNA.load,
            aes(x = Time,
                y = mean.viralLoad)) +
  geom_line(data = summary.vRNA.load[!is.na(summary.vRNA.load$mean.fragments),],
            aes(x = Time,
                y = (mean.fragments) * 5 * 10 ** 6),
            #color = "#1D7874",
            color = "gray20",
            #color = "yellow",
            #linewidth = 1.4,
            linetype = 2) +
  #scale_y_log10() +
  # Customize theme
  theme_classic() +
  theme(#axis.text =  element_blank(),
        #axis.title = element_blank()
  )


ggplot(vRNA.load,
       aes(x = Time,
           #y = Viral.load
           )
           ) +
  # geom_point(#size = 1
  # ) +
  geom_point(aes(y = mapped_viral_reads_percent * 5 * 10 ** 6),
             shape = 1,
             #fill = "#1D7874",
             #color = "#5F0F40", 
             color = "gray20",
             #color = "yellow",
             #size = 1
  ) +
  # mean values (3 rep) as lines
  # geom_line(data = summary.vRNA.load,
  #           aes(x = Time,
  #               y = mean.viralLoad)) +
  geom_line(data = summary.vRNA.load[!is.na(summary.vRNA.load$mean.fragments),],
            aes(x = Time,
                y = (mean.fragments) * 5 * 10 ** 6),
            #color = "#1D7874",
            color = "gray20",
            #color = "yellow",
            #linewidth = 1.4,
            linetype = 2) +
  xlab("Time (hpi)") +
  ylab("Viral load") +
  theme_classic() +
  theme(#axis.text =  element_blank(),
        axis.title = element_blank()
  )

# C) %infected worms
smF %>%
  pivot_longer(.,
               cols = c("RNA1", "RNA2", "lumen"),
               names_to = "signal",
               values_to = "% infected") %>%
  ggplot(.,
         aes(x = Time,
             y = `% infected`)
         ) +
  # geom_point(aes(color = signal),
  #            show.legend = F) +
  geom_jitter(aes(color = signal),
             show.legend = F,
             width = 0.10) +
  # mean values
  geom_line(data = means.smF,
            aes(x = Time,
                y = mean.RNA1),
            color = "#1CDFDC") +
  geom_line(data = means.smF,
            aes(x = Time,
                y = mean.RNA2),
            color = "#DB00DD") +
  geom_line(data = means.smF,
            aes(x = Time,
                y = mean.lumen)) +
  scale_color_manual(values = c("RNA1" = "#1CDFDC",
                                "RNA2" = "#DB00DD",
                                "lumen" = "black")) +
  xlim(c(0, 44)) +
  theme_classic() +
  theme(axis.text =  element_blank(),
        axis.title = element_blank()
  )
  

# E) ratio RNA2/RNA1 previous relativization to each fragment size
ggplot(infoSamples.complete[infoSamples.complete$Group == "infection",],
       aes(x = as.numeric(as.character(Time)),
           y = ratioRNA2.RNA1.sizeRelat)) +
  geom_point() +
  geom_line(data = mean.RNA.ratio,
            aes(x = Time,
                y = mean.ratio)) +
  
  # geom_line(data = mean.RNA.ratio,
  #           aes(x = Time,
  #               y = median.ratio),
  #           color = "red") +
  geom_hline(yintercept = 1,
             linetype = 2) +
  xlab("Time (hpi)") +
  ylab("RNA2 / RNA1") +
  theme_classic() +
  theme(#axis.text =  element_blank(),
        axis.title = element_blank()
  )
  ### probar a separa por proteínas?

## Reconfirm ratio RNA2/RNA1 are relative to the fragment length
infoSamples.complete <- infoSamples.complete %>%
  mutate(ratioRNA2.RNA1.sizeRelat.new = (RNA2_OrV_numReads/length(orv$ORV_RNA2_genome_our_mut)) / (RNA1_OrV_numReads/length(orv$ORV_RNA1_genome_our_mut)))

### MEAN AND SD RATIO RNA2/RNA1
mean(infoSamples.complete$ratioRNA2.RNA1.sizeRelat, na.rm = T)
sd(infoSamples.complete$ratioRNA2.RNA1.sizeRelat, na.rm = T)
  
vRNA.load$Time
max(vRNA.load$mapped_viral_reads_percent, na.rm = T)
infoSamples.complete$RNA2_OrV_numReads



ggplot(infoSamples.complete %>%
         filter(Group == "infection"),
       aes( x = RNA1_OrV_numReads /length(orv$ORV_RNA1_genome_our_mut),
            y = RNA2_OrV_numReads /length(orv$ORV_RNA2_genome_our_mut),
            group = Replica,
            color = Time)) +
  geom_path() +
  stat_summary(fun = mean,
             geom = "point") +
  #scale_y_log10() +
  #scale_x_log10() +
  #geom_point() +
  theme_classic()

### CORRELATION VIRAL RNA BY RT-QPCR & RNA-SEQ
summary_vRNA <- vRNA.load[complete.cases(vRNA.load),] %>%
  group_by(Time) %>%
  summarise(mean_vl.RTqPCR = mean(Viral.load),
            sd_vl.RTqPCT = sd(Viral.load),
            mean_vl.RNAseq = mean(mapped_viral_reads_percent),
            sd_vl.RNAseq = sd(mapped_viral_reads_percent))

# Means doesn't follow a normal distribution so --> Spearman's correlation
cor.test(summary_vRNA$mean_vl.RTqPCR,
         summary_vRNA$mean_vl.RNAseq,
         method = "spearman")
# Spearman's rank correlation rho
# 
# data:  summary_vRNA$mean_vl.RTqPCR and summary_vRNA$mean_vl.RNAseq
# S = 8, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.972028 

