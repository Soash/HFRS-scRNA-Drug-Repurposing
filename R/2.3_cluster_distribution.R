library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)

# ==================================================================
# Environment Setup
# ==================================================================
rm(list = ls())
gc(); gc()

setwd("~/HTNV")

# Create output directory
result_dir <- "R/result/2.3_distribution"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

### ==================================================================
### PHASE 2.3: CALCULATE PERCENTAGES AND PLOT DISTRIBUTIONS
### ==================================================================
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading Annotated Object...\n")
combined <- qs_read("R/result/2.2_annotation/2.2_annotation.qs2")

# Extract metadata
meta <- combined@meta.data

cat(format(Sys.time(), "[%H:%M:%S]"), "Step 2: Calculating Cell Counts and Percentages...\n")

# Calculate absolute counts per cell type and condition
count_data <- meta %>%
  group_by(cell_type, condition) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate Percentages: "What percentage of the Cell Type is Normal vs Disease?"
pct_within_cluster <- count_data %>%
  group_by(cell_type) %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

# --- EXPORT TO CSV ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 3: Exporting CSV...\n")

csv_output <- pct_within_cluster %>%
  select(cell_type, condition, percentage) %>%
  pivot_wider(names_from = condition, values_from = percentage, values_fill = 0) %>%
  # Ensure columns exist even if one condition is missing in some clusters
  mutate(
    Normal_Percent = if("Normal" %in% names(.)) Normal else 0,
    Disease_Percent = if("Disease" %in% names(.)) Disease else 0
  ) %>%
  select(Cluster_Cell_Type = cell_type, Healthy_Percent = Normal_Percent, Disease_Percent = Disease_Percent) %>%
  # Round to 2 decimal places for cleaner viewing
  mutate(across(where(is.numeric), ~round(., 2)))

write.csv(csv_output, file.path(result_dir, "Cluster_Healthy_Disease_Percentages.csv"), row.names = FALSE)

# --- GENERATE PLOTS ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 4: Generating Plots...\n")

# Plot 1: Stacked Bar Plot (Proportion of Disease vs Normal within each Cell Type)
# This perfectly visualizes your CSV data.
plot1 <- ggplot(meta, aes(x = cell_type, fill = condition)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Disease" = "#E41A1C")) +
  labs(
    title = "Proportion of Healthy vs Disease Cells per Cluster",
    x = "Cell Type / Cluster",
    y = "Percentage",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Plot 2: Total Cell Counts (To see if a 100% disease cluster is just 5 cells or 5000 cells)
plot2 <- ggplot(meta, aes(x = cell_type, fill = condition)) +
  geom_bar(position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Disease" = "#E41A1C")) +
  labs(
    title = "Absolute Cell Counts per Cluster",
    x = "Cell Type / Cluster",
    y = "Number of Cells",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save the plots
pdf(file.path(result_dir, "Cluster_Distribution_Plots.pdf"), width = 12, height = 7)
print(plot1)
print(plot2)
dev.off()

cat(format(Sys.time(), "[%H:%M:%S]"), "COMPLETE! Check the 'R/result/2.3_distribution' folder.\n")