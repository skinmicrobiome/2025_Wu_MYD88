library(ggplot2)
library(dplyr)
library(ggsci)
library(ggforce)
library(vegan)
library(svglite)

abun <- read.csv('fig5bd.csv')
abun$label <- factor(abun$label, levels = c(
  'Actinomycetota', 'Corynebacterium',
  'Bacillota', 'Staphylococcus', 'Streptococcus',
  'Pseudomonadota',
  'Bacteroidota', 'Bacteroides', 'Alistipes',
  'Verrucomicrobiota', 'Akkermansia',
  'Bacteria (Other)'
))
abun$treatment_group <- factor(abun$treatment_group, levels = c('ft/ft', 'veh', 'neo'))

custom_palette <- c(
  "Bacillota"         = "#c6dbef",
  "Staphylococcus"    = "#5c82b2",
  "Streptococcus"     = "#008ECC",
  "Pseudomonadota"    = "#FF99CC",
  "Actinomycetota"    = "#fee6ce",
  "Corynebacterium"   = "#e6550d",
  "Bacteroidota"      = "#c7e9c0",
  "Bacteroides"       = "#31a354",
  "Alistipes"         = "#66c2a4",
  "Verrucomicrobiota" = "#81C7D4",
  "Akkermansia"       = "#a1d8e6",
  "Bacteria (Other)"  = "#d9d9d9"
)

bar_theme <- theme(
  panel.spacing   = unit(0, "lines"),
  text            = element_text(size = 14),
  axis.text       = element_text(size = 14),
  axis.title.x    = element_text(size = 14),
  axis.title.y    = element_text(size = 18),
  legend.text     = element_text(size = 14),
  legend.title    = element_blank(),
  legend.key.size = unit(1.2, "cm"),
  strip.text.x    = element_text(size = 16),
  strip.text.y    = element_text(size = 16),
  axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5)
)

# Fig5b: Day 0
d0 <- abun %>% filter(timepoint == 'Day 0')

plot1 <- ggplot(d0, aes(x = mouseID, y = rel_abun, fill = label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  ylab('Relative abundance %') +
  xlab('') +
  facet_grid(timepoint ~ index, scales = "free", space = 'free') +
  bar_theme

ggsave("fig5b.svg", plot = plot1, width = 12, height = 6, units = "in")

# Fig5d: Day 7
d7 <- abun %>% filter(timepoint == 'Day 7')

plot2 <- ggplot(d7, aes(x = mouseID, y = rel_abun, fill = label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  ylab('Relative abundance %') +
  xlab('') +
  facet_grid(timepoint ~ treatment_group, scales = "free", space = 'free') +
  bar_theme

ggsave("fig5d.svg", plot = plot2, width = 12, height = 6, units = "in")

# FigS2
my_data     <- read.csv("figs2.csv", header = TRUE, row.names = 1)
data_matrix <- my_data[, 1:(ncol(my_data) - 3)]
bc_dist     <- vegdist(data_matrix, method = "bray")
pcoa        <- cmdscale(bc_dist, k = 2, eig = TRUE)
eig         <- pcoa$eig
perc1       <- eig[1] / sum(eig) * 100
perc2       <- eig[2] / sum(eig) * 100

pcoa_df <- data.frame(
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  group = factor(my_data$group, levels = c("ft/ft", "neo.Day 0", "veh.Day 7", "neo.Day 7"))
)

group_colors <- c(
  "ft/ft"     = "#F8766D",
  "neo.Day 0" = "#7CAE00",
  "neo.Day 7" = "#00BFC4",
  "veh.Day 7" = "#C77CFF"
)

p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = group, color = group)) +
  stat_ellipse(aes(group = group),
               alpha = 1, linetype = "dashed", linewidth = 1.5) +
  geom_point(shape = 21, color = "black", size = 3, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", perc1),
    y = sprintf("PCoA2 (%.1f%%)", perc2)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title    = element_blank(),
    text            = element_text(size = 16),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 16),
    panel.border    = element_blank(),
    axis.line       = element_blank()
  ) +
  guides(fill  = guide_legend(nrow = 1, byrow = TRUE),
         color = guide_legend(nrow = 1, byrow = TRUE))

ggsave("figs2.svg", plot = p, width = 7, height = 7, dpi = 300)