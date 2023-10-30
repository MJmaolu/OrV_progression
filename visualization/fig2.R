################################################################################
# LOCALIZATION OF vRNA IN INFECTED WORMS BY TIME
#
# MJ Olmo-Uceda
# 2023/08/01
################################################################################
# fig2
################################################################################
localization.infection

localization.infection$`Cell&Lumen` <- localization.infection$CellPartialLumen + localization.infection$CellTotalLumen

# Grouping in only 3 categories: Cell, Cell&lumen and Lumen
localization.infection$Lumen <- localization.infection$PartialLumen + localization.infection$TotalLumen
localization.infection.long <- localization.infection %>%
  dplyr::select(c("time", "Cell", "Cell&Lumen", "Lumen")) %>%
  pivot_longer(.,
               cols = c("Cell", "Cell&Lumen", "Lumen"),
               names_to = "localization vRNA",
               values_to = "percentage infected worms")

localization.infection.long$`localization vRNA` <- factor(localization.infection.long$`localization vRNA`,
                                                          levels = c("Cell", "Cell&Lumen", "Lumen"))


ggplot(localization.infection.long, 
       aes(x = factor(time), 
           y = `percentage infected worms`,
           fill = `localization vRNA`)) +
  geom_bar(stat = "identity",
           #position = "dodge",
           width = 1,
           color = "black") +
  labs(#title = "Composition Plot",
       x = "Time",
       y = "Localization (%)") +
  theme_minimal() +
  scale_fill_manual(values = c("Cell" = "#A3D2C6",#"#E7C633",
                               "Cell&Lumen" = "#FFFDB5",#DD625E",
                               "Lumen" = "#BBB8D7"#3E474B"
                               )) +
  theme(axis.text = element_text(color = "black",
                                 size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "black",
                                          linewidth = 0.2),
        panel.grid.minor.y = element_blank(),
        axis.title = element_blank()
        )
################################################################################
# fig2.R
# MJ Olmo-Uceda
# 2023/04/21
################################################################################
str(area.inf)

ggplot(area.inf, 
       aes(y = idGus)) +
  geom_segment(aes(x = X,
                   y = Y,
                   xend = xend,
                   yend = yend,
                   color = "#000000",
  ),
  linewidth = 2.2,
  lineend = "round"
  ) +
  geom_segment(aes(x = X,
                   y = Y,
                   xend = xend,
                   yend = yend,
                   color = as.factor(time),
                   ),
               linewidth = 1.5,
               lineend = "round"
               ) +
  xlab("Relative worm length") +
  ylab("") +
  
  xlim(c(4,96)) +
  facet_wrap(~ time,
             nrow = 2,
             scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("10" = "#D07DAB",#FBBFE0",
                                "24" = "#B7DBA9" #9DD893"
                                  )) +
  theme(axis.text = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text = element_blank()
        )
str(area.inf)


# B) Number infected cells per hpi
str(num.cells.hour)
num.cells.hour$hours <- factor(num.cells.hour$hours,
                               levels = unique(num.cells.hour$hours))

ggplot(num.cells.hour[num.cells.hour$cells != 0,],
       aes(y = fct_rev(hours), # up to down
           x = as.numeric(cells),
           fill = hours)) +
  geom_density_ridges2(aes(fill = hours),
                       stat = "binline",
                       binwidth = 1,
                       scale = 0.95) +
  theme_classic() +
  scale_fill_manual(values = c(rep("#D07DAB",5), 
                               rep("#E9F4A3", 1),
                               rep("#B7DBA9", 2),
                               rep("#536A75", 2)) 
  ) +
  scale_y_discrete(expand = c(0,0)) +
  # scale_x_discrete(breaks = seq(1, 
  #                               max(n_cells_mutants$Cells), 
  #                               by = 1)) +  # Set tick breaks
  scale_x_continuous(
    breaks = c(0:5),
    limits = c(0, 5),
    expand = c(0, 0),
    name = "Number of infected cells"
  ) +
  xlab("Number of infected cells") +
  ylab("") +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(15, 15, 15, 15, "pt")
  )

num.cells.hour[num.cells.hour$cells != 0,] %>%
  group_by(hours, cells) %>%
  dplyr::count(cells) %>%
  ggplot(.,
         aes(y = fct_rev(hours),
             x = as.numeric(cells),
             color = hours,
             size = n*10)) +
  geom_point() +
  scale_color_manual(values = c(rep("#E77519",5), 
                               rep("#64C28D", 1),
                               rep("#AF00DE", 2),
                               rep("black", 2))) +
  theme_classic()
    
  

