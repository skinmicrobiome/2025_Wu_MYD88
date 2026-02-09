library(vegan)
library(ggplot2)
library(dplyr)
library(svglite)

my_data <- read.csv("figs5.csv", header = TRUE, row.names = 1)

data_matrix <- my_data[, 1:(ncol(my_data) - 3)]
# Bray-Curtis distance
bc_dist <- vegdist(data_matrix, method = "bray")

# PCoA
pcoa  <- cmdscale(bc_dist, k = 2, eig = TRUE)
eig   <- pcoa$eig
perc1 <- eig[1] / sum(eig) * 100
perc2 <- eig[2] / sum(eig) * 100
# Plotting dataframe
pcoa_df <- data.frame(
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  group = factor(my_data$group, levels = c("ft/ft", "ft/ftXMyD88-/-", "ft/ftXMyD88-/- veh", "ft/ftXMyD88-/- neo"))
)

group_colors <- c(
  "ft/ft" = "#F8766D",
  "ft/ftXMyD88-/-" = "#7CAE00", 
  "ft/ftXMyD88-/- neo" = "#00BFC4",
  "ft/ftXMyD88-/- veh" = "#C77CFF"
)
p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = group, color = group)) +
  stat_ellipse(aes(group = group),
               alpha = 1, linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(fill = group), color = "black", size = 3, alpha = 0.7, shape = 21) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", perc1),
    y = sprintf("PCoA2 (%.1f%%)", perc2)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    panel.border = element_blank(),
    axis.line = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         color = guide_legend(nrow = 1, byrow = TRUE))
ggsave("figS5.svg", plot = p, width = 7, height = 7, dpi = 300)
#print(p)