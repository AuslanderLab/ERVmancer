#This is an Rscript to plot the output of the results for the figure 4.

library(ggplot2)
library(ggpubr)
library(ggrepel)



#First, look at all of the combined ones------------------
all_combined <- read.csv("data_for_R/combined_all_mwu_sig_sum_counts_for_R.csv")

# Set p53_status factor order
all_combined$p53_status <- factor(all_combined$p53_status, levels = c("wt", "del", "mut"))

# Define comparisons
comparisons <- list(
  c("wt", "del"),
  c("wt", "mut"),
  c("del", "mut")
)

plot <- ggplot(all_combined, aes(x = p53_status, y = normalized_counts, fill = p53_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(color = "black", width = 0.2, size = 2, alpha = 0.6) +  # Black points
  stat_compare_means(
    comparisons = comparisons,
    label = "p.format",
    method = "wilcox.test",
    label.y = c(180, 200, 235)
  ) +
  coord_cartesian(ylim = c(0, 265)) +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +  # Box colors
  labs(
    x = NULL,
    y = "Reads Per Million"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    legend.position = "none",
    axis.ticks.y = element_line(color = "black"),
    panel.background = element_blank()
  )

# Save the plot as PDF
ggsave("R_figures/all_combined_plot.pdf", plot = plot, width = 8, height = 6)

# Save the plot as JPG
ggsave("R_figures/all_combined_plot.jpg", plot = plot, width = 8, height = 6, dpi = 300)


#second, look at just ERV_106-------------------------------

erv_106 <- read.csv('data_for_R/ERV_1069516_normcounts_for_R.csv')

# Set factor order
erv_106$p53_status <- factor(erv_106$p53_status, levels = c("wt", "del", "mut"))

# Define comparisons
comparisons <- list(
  c("wt", "del"),
  c("wt", "mut"),
  c("del", "mut")
)

plot <- ggplot(erv_106, aes(x = p53_status, y = ERV_1069516_HERVH_int_LTR7Y, fill = p53_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Hide outliers so we can control point appearance manually
  geom_jitter(color = "black", width = 0.2, size = 2, alpha = 0.6) +  # Add visible black points
  stat_compare_means(
    comparisons = comparisons,
    label = "p.format",
    method = "wilcox.test",
    label.y = c(1.6, 1.8, 2.0)
  ) +
  coord_cartesian(ylim = c(0, 2.1)) +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  labs(
    x = NULL,
    y = "Reads Per Million"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    legend.position = "none",
    axis.ticks.y = element_line(color = "black"),
    panel.background = element_blank()
  )
# Save the plot as PDF
ggsave("R_figures/erv_106_plot.pdf", plot = plot, width = 8, height = 6)

# Save the plot as JPG
ggsave("R_figures/erv_106_plot.jpg", plot = plot, width = 8, height = 6, dpi = 300)



#--------------------------------------------
#Third, look at ERV_174 (actually clade_2844, from which we extracted ERV_174)-------

erv_174 <- read.csv('data_for_R/Clade_2844_normcounts_for_R.csv')

# Set p53_status as a factor with the desired order
erv_174$p53_status <- factor(erv_174$p53_status, levels = c("wt", "del", "mut"))

# Define comparisons
comparisons <- list(
  c("wt", "del"),
  c("wt", "mut"),
  c("del", "mut")
)

plot <- ggplot(erv_174, aes(x = p53_status, y = Clade_2844, fill = p53_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Hide outliers
  geom_jitter(color = "black", width = 0.2, size = 2, alpha = 0.6) +  # Add black points
  
  # Display exact p-values
  stat_compare_means(
    comparisons = comparisons,
    label = "p.format",
    method = "wilcox.test",
    label.y = c(45, 55, 65)
  ) +
  
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  labs(
    x = NULL,
    y = "Normalized Counts"
  ) +
  coord_cartesian(ylim = c(0, 75)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 35),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    axis.ticks.y = element_line(color = "black"),
    panel.background = element_blank()
  )

# Save the plot as PDF
ggsave("R_figures/clade_2844_plot.pdf", plot = plot, width = 8, height = 6)

# Save the plot as JPG
ggsave("R_figures/clade_2844_plot.jpg", plot = plot, width = 8, height = 6, dpi = 300)

