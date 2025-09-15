#This is an R file to try and make the figures for the paper.

library(dplyr)
library('ggplot2')
library(tidyr)
library(ggExtra)

#figure 2A---------------------------------------------------

corr_df = read.csv('data_for_R/fig2A_correlation_df.csv',header = TRUE, sep = ',')

corr_df$X <- reorder(corr_df$X, -corr_df$statistic)

fig2a <- ggplot(corr_df, aes(x = X, y = statistic, fill = X)) +
  geom_bar(stat = "identity") +
  labs(
    x = NULL, 
    y = "Long Read SpearmanR Correlation"
  ) +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  scale_y_continuous(
    breaks = seq(0, 0.8, by = 0.1),  # Set y-axis increments to 0.1
    limits = c(0, 0.8)               # Set the y-axis range to 0-0.8
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),  # Rotate x-axis labels without changing size or boldness
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    axis.ticks.y = element_line(color = "black")
  )
fig2a
# Export the plot to a PDF
ggsave("R_figures/figure2A.pdf", plot = fig2a, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure2A.jpg", plot = fig2a, device = "jpg", width = 4, height = 6)


#Figure 2B,C---------------------------------------------------
bootstrap_data = read.csv('data_for_R/fig2B_melted_bootstrap_data_80.csv', header = TRUE, sep = ',')

# Ensure Columns in bootstrap_data follows the same order as X in corr_df
bootstrap_data$Columns <- factor(bootstrap_data$Columns, levels = levels(corr_df$X))

# Replot with the correct order
fig2b <- ggplot(bootstrap_data, aes(x = Columns, y = Values, fill = Columns)) +
  geom_boxplot(size = 0.2, outlier.size = 0.9) +  # Make the box outlines thinner
  labs(x = NULL, y = "Long Read SpearmanR Correlation") +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),  # Rotate x-axis labels without changing size or boldness
    legend.position = "none"
  )

print(fig2b)

ggsave("R_figures/figure2B_80.pdf", plot = fig2b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure2B_80.jpg", plot = fig2b, device = "jpg", width = 4, height = 6)

#Do the same thing for the 50 percent downsampling

bootstrap_data = read.csv('data_for_R/fig2C_melted_bootstrap_data_50.csv', header = TRUE, sep = ',')

# Ensure Columns in bootstrap_data follows the same order as X in corr_df
bootstrap_data$Columns <- factor(bootstrap_data$Columns, levels = levels(corr_df$X))


fig2c <- ggplot(bootstrap_data, aes(x = Columns, y = Values, fill = Columns)) +
  geom_boxplot(size = 0.2, outlier.size = 0.9) +  # Make the box outlines thinner
  labs(x = NULL, y = "Long Read SpearmanR Correlation") +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),  # Rotate x-axis labels without changing size or boldness
    legend.position = "none"
  )

print(fig2c)

ggsave("R_figures/figure2C_50.pdf", plot = fig2c, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure2C_50.jpg", plot = fig2c, device = "jpg", width = 4, height = 6)



#Figure 2D------------------------------------
single_sample_corr = read.csv('data_for_R/fig2D_correlation_df_by_sample.csv', header = TRUE,sep = ',')
single_sample_corr$Variable <- reorder(single_sample_corr$Variable, single_sample_corr$Value, FUN = mean, decreasing = TRUE)
# Add group (H, A, S) column
single_sample_corr <- single_sample_corr %>%
  mutate(group = factor(substr(sample, 1, 1), levels = c("H", "A", "S")))

# Order sample levels
single_sample_corr$sample <- factor(
  single_sample_corr$sample,
  levels = c(
    grep("^H", unique(single_sample_corr$sample), value = TRUE),
    grep("^A", unique(single_sample_corr$sample), value = TRUE),
    grep("^S", unique(single_sample_corr$sample), value = TRUE)
  )
)

# Calculate mean and standard error
mean_se_data <- single_sample_corr %>%
  group_by(Variable) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define fill colors for bars (by Variable)
bar_colors <- c("#457B9D", "#F4A261", "#E9C46A")  # original bar colors

# Define point colors for group (H, A, S)
group_colors <- c(
  H = "#B22222",  # darker red (firebrick)
  A = "#006400",  # very dark green
  S = "#4B0082"   # deep indigo (contrasts well with yellow/orange)
)

# Count how many samples per Variable
single_sample_corr <- single_sample_corr %>%
  group_by(Variable) %>%
  arrange(group, .by_group = TRUE) %>%
  mutate(n_in_group = n(),
         row_number = row_number(),
         x_jittered = as.numeric(Variable) + 
           (row_number - (n_in_group + 1)/2) * 0.05)  # adjust 0.05 for spacing


# Plot
fig2d <- ggplot(mean_se_data, aes(x = Variable, y = mean_value, fill = Variable)) +  
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +
  geom_point(
    data = single_sample_corr,
    aes(x = x_jittered, y = Value, color = group),
    size = 2, alpha = 1
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_color_manual(
    name = "Sample Group",
    values = group_colors,
    labels = c("Healthy Control", "Active MS", "Stable MS")
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 8, by = 0.1)) +
  labs(y = "Long Read Spearman Correlation") +
  guides(fill = "none") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    legend.position = "top",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    axis.ticks.y = element_line(color = "black")
  )

fig2d

ggsave("R_figures/figure2d.pdf", plot = fig2d, device = "pdf", width = 6, height = 6)

ggsave("R_figures/figure2d.jpg", plot = fig2d, device = "jpg", width = 6, height = 6)


#Figure 2E-G------------------------------------------------
#These figures are the correlation scatterplots

df1 <- read.csv('data_for_R/fig2E_ervmancer_against_minimap_corrplot.csv')

# Define custom colors
custom_colors <- c(
  "#E69F00", "#D55E00", "#F0E442",  # Three shades of orange for H
  "#56B4E9", "#0072B2", "#6699CC",  # Three shades of blue for A
  "#009E73", "#117733", "#44AA99", "#88CCAA"  # Four shades of green for S
)


# Reorder the 'sample' factor based on groups
df1$sample <- factor(
  df1$sample,
  levels = c(
    grep("^H", unique(df1$sample), value = TRUE),
    grep("^A", unique(df1$sample), value = TRUE),
    grep("^S", unique(df1$sample), value = TRUE)
  )
)

# Base scatter plot
fig2e <- ggplot(df1, aes(x = cpm1, y = cpm2, color = sample)) +
  geom_point(size = 0.5, alpha = 0.4) +  # Smaller points
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = NULL,
    x = "ERVmancer -log",
    y = "minimap2 -log",
    color = "Cell Line"  # Legend title
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1),  # Center title
    legend.position = "right",  # Move legend to the right
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 14)  # Increase axis label size
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1))  # Larger dots in the legend
  )

# Add marginal density plots
fig2e_with_marginals <- ggMarginal(
  fig2e,
  type = "density",    # Use density plots
  margins = "both",    # Add distributions to both x and y axes
  size = 5,            # Size of the marginal plots
  fill = "gray80",     # Fill under the density lines
  color = "black"      # Outline color for the density plots
)

fig2e_with_marginals

ggsave("R_figures/figure2e.pdf", plot = fig2e_with_marginals, device = "pdf", width = 6, height = 6)
#ggsave("R_figures/figure2d.jpg", plot = fig2d_with_marginals, device = "jpg", width = 6, height = 6)
ggsave("R_figures/figure2e.jpg", 
       plot = fig2e_with_marginals, 
       device = "jpg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  # Use dpi instead of res



#Figure 3F-------------------------------------------------
df1 <- read.csv('data_for_R/fig2F_ervmap_against_minimap_corrplot.csv')

# Define custom colors
custom_colors <- c(
  "#E69F00", "#D55E00", "#F0E442",  # Three shades of orange for H
  "#56B4E9", "#0072B2", "#6699CC",  # Three shades of blue for A
  "#009E73", "#117733", "#44AA99", "#88CCAA"  # Four shades of green for S
)

# Reorder the 'sample' factor based on groups
df1$sample <- factor(
  df1$sample,
  levels = c(
    grep("^H", unique(df1$sample), value = TRUE),
    grep("^A", unique(df1$sample), value = TRUE),
    grep("^S", unique(df1$sample), value = TRUE)
  )
)

# Base scatter plot
fig2f <- ggplot(df1, aes(x = cpm1, y = cpm2, color = sample)) +
  geom_point(size = 0.5, alpha = 0.4) +  # Smaller points
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = NULL,
    x = "ERVmap -log",
    y = "minimap2 -log",
    color = "Cell Line"  # Legend title
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1),  # Center title
    legend.position = "right",  # Move legend to the right
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 14)  # Increase axis label size
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1))  # Larger dots in the legend
  )

# Add marginal density plots
fig2f_with_marginals <- ggMarginal(
  fig2f,
  type = "density",    # Use density plots
  margins = "both",    # Add distributions to both x and y axes
  size = 5,            # Size of the marginal plots
  fill = "gray80",     # Fill under the density lines
  color = "black"      # Outline color for the density plots
)

fig2f_with_marginals

ggsave("R_figures/figure2F.pdf", plot = fig2f_with_marginals, device = "pdf", width = 6, height = 6)
#ggsave("R_figures/figure2F.jpg", plot = fig2d_with_marginals, device = "jpg", width = 6, height = 6)
ggsave("R_figures/figure2F.jpg", 
       plot = fig2f_with_marginals, 
       device = "jpg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  # Use dpi instead of res



#Figure 3G-------------------------------------------------
df1 <- read.csv('data_for_R/fig2G_telescope_against_minimap_corrplot.csv')

# Define custom colors
custom_colors <- c(
  "#E69F00", "#D55E00", "#F0E442",  # Three shades of orange for H
  "#56B4E9", "#0072B2", "#6699CC",  # Three shades of blue for A
  "#009E73", "#117733", "#44AA99", "#88CCAA"  # Four shades of green for S
)

# Reorder the 'sample' factor based on groups
df1$sample <- factor(
  df1$sample,
  levels = c(
    grep("^H", unique(df1$sample), value = TRUE),
    grep("^A", unique(df1$sample), value = TRUE),
    grep("^S", unique(df1$sample), value = TRUE)
  )
)

# Base scatter plot
fig2g <- ggplot(df1, aes(x = cpm1, y = cpm2, color = sample)) +
  geom_point(size = 0.5, alpha = 0.4) +  # Smaller points
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(
    title = NULL,
    x = "TELEscope -log",
    y = "minimap2 -log",
    color = "Cell Line"  # Legend title
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 1),  # Center title
    legend.position = "right",  # Move legend to the right
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 14)  # Increase axis label size
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1))  # Larger dots in the legend
  )

# Add marginal density plots
fig2g_with_marginals <- ggMarginal(
  fig2g,
  type = "density",    # Use density plots
  margins = "both",    # Add distributions to both x and y axes
  size = 5,            # Size of the marginal plots
  fill = "gray80",     # Fill under the density lines
  color = "black"      # Outline color for the density plots
)

fig2g_with_marginals

ggsave("R_figures/figure2G.pdf", plot = fig2g_with_marginals, device = "pdf", width = 6, height = 6)
#ggsave("R_figures/figure2G.jpg", plot = fig2d_with_marginals, device = "jpg", width = 6, height = 6)
ggsave("R_figures/figure2G.jpg", 
       plot = fig2g_with_marginals, 
       device = "jpg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  # Use dpi instead of res

