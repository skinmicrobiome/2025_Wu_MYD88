# libraries
library(ggplot2)
library(dplyr)
library(ggsci)
library(ggforce)

abun <- read.csv('fig4bd.csv')

# taxonomic order and treatment group order
abun$label <- factor(abun$label, levels = c(
  'Actinobacteriota', 'Corynebacterium',
  'Firmicutes', 'Staphylococcus', 'Streptococcus',
  'Proteobacteria',
  'Bacteroidota', 'Bacteroides', 'Alistipes',
  'Verrucomicrobiota', 'Akkermansia',
  'Bacteria (Other)'
))
abun$treatment_group <- factor(abun$treatment_group, levels = c('ft/ft', 'veh', 'neo'))

# custom color palette
custom_palette <- c(
  "Firmicutes" = "#c6dbef",
  "Staphylococcus" = "#5c82b2",
  "Streptococcus" = "#008ECC",
  "Proteobacteria" = "#FF99CC",
  "Actinobacteriota" = "#fee6ce",
  "Corynebacterium" = "#e6550d",
  "Bacteroidota" = "#c7e9c0",
  "Bacteroides" = "#31a354",
  "Alistipes" = "#66c2a4",
  "Verrucomicrobiota" = "#81C7D4",
  "Akkermansia" = "#a1d8e6",
  "Bacteria (Other)" = "#d9d9d9"
)

# Day 0
d0 <- abun %>% filter(timepoint == 'Day 0')

plot1 <- ggplot(d0, aes(x = mouseID, y = rel_abun, fill = label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  ylab('Relative abundance %') +
  xlab('') +
  facet_grid(timepoint ~ index, scales = "free", space = 'free') +
  theme(
    panel.spacing = unit(0, "lines"), 
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.key.size = unit(1.2, "cm"),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("fig4b.svg", plot = plot1, width = 12, height = 6, units = "in")

# Day 7 
d7 <- abun %>% filter(timepoint == 'Day 7')

plot2 <- ggplot(d7, aes(x = mouseID, y = rel_abun, fill = label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  ylab('Relative abundance %') +
  xlab('') +
  facet_grid(timepoint ~ treatment_group, scales = "free", space = 'free') +
  theme(
    panel.spacing = unit(0, "lines"), 
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.key.size = unit(1.2, "cm"),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("fig4d.svg", plot = plot2, width = 12, height = 6, units = "in")
