#This is an R file to make the figure 3 for the paper.

library(dplyr)
library('ggplot2')
library(tidyr)


#Figure 3A-------------------------------------------------------

#This is correlations with simulated data
corr_df = read.csv('data_for_R/correlation_df_10_samples.csv',header = TRUE, sep = ',')
corr_df
# Reshape data from wide to long format
corr_long <- corr_df %>%
  pivot_longer(cols = -X, names_to = "Metric", values_to = "Value")

# Compute mean, standard deviation, and standard error for each metric
error_data <- corr_long %>%
  group_by(Metric) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    sem_value = sd_value / sqrt(n()),  # Standard Error of the Mean (SEM)
    ci_95 = 1.96 * sem_value           # 95% Confidence Interval
  )

# Merge with original data
corr_long <- corr_long %>%
  left_join(error_data, by = "Metric")

# Create the plot with error bars
fig2a <- ggplot(corr_long, aes(x = Metric, y = mean_value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  
  geom_errorbar(aes(ymin = mean_value - ci_95, ymax = mean_value + ci_95), width = 0.2) +  
  labs(
    x = NULL, 
    y = "Ground Truth SpearmanR Correlation"
  ) +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),  # Keep 0.1 increments
    limits = c(0, 1)               # Force the y-axis to extend to 1
  ) +  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),  
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    axis.ticks.y = element_line(color = "black")
  )

# Print the plot
print(fig2a)

# Export the plot to a PDF
ggsave("R_figures/figure3A.pdf", plot = fig2a, device = "pdf", width = 4, height = 6)

ggsave("R_figures/figure3A.jpg", plot = fig2a, device = "jpg", width = 4, height = 6)

#Figure 3B---------------------------------------------------
#this is code to make figure 3B, which is looking at the False Positive Rate

input_data <- read.csv('data_for_R/fpr_df.csv', header = TRUE, sep = ',')
#input_data <- read.csv('data_for_analysis/fpr_cutoff_50_df_for_plotting.csv', header = TRUE, sep = ',')

# Transform data into long format
input_data_long <- input_data %>%
  pivot_longer(cols = -X, names_to = "Condition", values_to = "Value") %>%
  mutate(Condition = factor(Condition, levels = colnames(input_data)[-1]))  # Preserve original order

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

max_y <- max(summary_data$Mean + summary_data$SE)

fig3a <- ggplot(summary_data, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "False Positive Rate") +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("ERVmancer", "ERVmap", "TELEscope")) +
  coord_cartesian(ylim = c(0, max_y * 1.05)) +  # extend y-axis a little past max
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )
fig3a <- fig3a + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

print(fig3a)

ggsave("R_figures/figure3B.pdf", plot = fig3a, device = "pdf", width = 4, height = 6)

ggsave("R_figures/figure3B.jpg", plot = fig3a, device = "jpg", width = 4, height = 6)


#Figure 3C--------------------------------------
# same thing with Recall
input_data <- read.csv('data_for_R/recall_df.csv', header = TRUE, sep = ',')


# Transform data into long format
input_data_long <- input_data %>%
  pivot_longer(cols = -X, names_to = "Condition", values_to = "Value") %>%
  mutate(Condition = factor(Condition, levels = colnames(input_data)[-1]))  # Preserve original order

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Create bar plot with error bars
fig3b <- ggplot(summary_data, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "Recall") +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("ERVmancer", "ERVmap", "TELEscope")) +
  scale_y_continuous(limits = c(0, 1))+ #I can use this because none of the values are above 1
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

print(fig3b)

ggsave("R_figures/figure3C.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3C.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)

#Figure 3D---------------------------------------------------------
#This is the F1 Score, as recommended by Avi.
input_data <- read.csv('data_for_R/f1_df.csv', header = TRUE, sep = ',')


# Transform data into long format
input_data_long <- input_data %>%
  pivot_longer(cols = -X, names_to = "Condition", values_to = "Value") %>%
  mutate(Condition = factor(Condition, levels = colnames(input_data)[-1]))  # Preserve original order

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Create bar plot with error bars
fig3b <- ggplot(summary_data, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "F1 Score") +
  scale_fill_manual(values = c("#457B9D", "#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("ERVmancer", "ERVmap", "TELEscope")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

print(fig3b)

ggsave("R_figures/figure3D.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3D.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)


#Figure 3E Stuff-----------------------------------------
#This next code creates the many graphs for figure 3E
input_data_long <- read.csv('data_for_R/df_long_ervmancer_actual_1409.csv')

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Find max value
y_max <- max(summary_data$Mean + summary_data$SE, na.rm = TRUE)

fig3b <- ggplot(summary_data, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "Counts") +
  scale_fill_manual(values = c("#2A9D8F","#457B9D")) +
  scale_x_discrete(labels = c("Actual Counts", "ERVmancer")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  expand_limits(y = y_max + 5)

print(fig3b)

ggsave("R_figures/figure3E_clade_1409.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3E_clade_1409.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)

#Next, need for each of the ervs, for actual, ervmap, and TELEscope-----------------

input_data_long <- read.csv('data_for_R/df_long_4317323.csv')

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Find max value
y_max <- max(summary_data$Mean + summary_data$SE, na.rm = TRUE)

fig3b <- ggplot(summary_data, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "Counts") +
  scale_fill_manual(values = c("#2A9D8F","#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("Actual Counts", "ERVmap","TELEscope")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  expand_limits(y = y_max + 5)

print(fig3b)

ggsave("R_figures/figure3E_erv_4317323.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3E_erv_4317323.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)

#Finally, the last erv---------------------------------------------

input_data_long <- read.csv('data_for_R/df_long_4409086.csv')

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Find max value
y_max <- max(summary_data$Mean + summary_data$SE, na.rm = TRUE)

fig3b <- ggplot(summary_data, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "Counts") +
  scale_fill_manual(values = c("#2A9D8F","#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("Actual Counts", "ERVmap","TELEscope")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  ylim(0, 30) 

print(fig3b)


ggsave("R_figures/figure3E_erv_4409086.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3E_erv_4409086.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)

#Figure 3F--------------------------------------------------
#Do essentially the same thing for figure 3F

input_data_long <- read.csv('data_for_R/df_long_ERV_1481301.csv')

# Summarize data to get mean and standard error
summary_data <- input_data_long %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

# Find max value
y_max <- max(summary_data$Mean + summary_data$SE, na.rm = TRUE)

fig3b <- ggplot(summary_data, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(x = NULL, y = "Counts") +
  scale_fill_manual(values = c("#2A9D8F","#457B9D","#F4A261", "#E9C46A")) +
  scale_x_discrete(labels = c("Actual Counts", "ERVmancer","ERVmap","TELEscope")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  expand_limits(y = y_max + 5)

print(fig3b)

ggsave("R_figures/figure3F_erv_1481301.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
#ggsave("R_figures/figure3B_50_cutoff.pdf", plot = fig3b, device = "pdf", width = 4, height = 6)
ggsave("R_figures/figure3F_erv_1481301.jpg", plot = fig3b, device = "jpg", width = 4, height = 6)





