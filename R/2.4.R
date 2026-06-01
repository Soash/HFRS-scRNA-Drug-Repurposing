library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)

rm(list = ls())
gc(); gc()

setwd("~/HTNV")

result_dir <- "R/result/2.3_distribution"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

cat(format(Sys.time(), "[%H:%M:%S]"), "Loading object...\n")

combined <- qs_read("R/result/2.2_annotation/2.2_annotation.qs2")
meta <- combined@meta.data

# ==========================================================
# Count cells per cell type and condition
# ==========================================================

count_data <- meta %>%
  group_by(cell_type, condition) %>%
  summarise(cell_count = n(), .groups = "drop")

# ==========================================================
# Total cells per condition
# ==========================================================

total_cells <- meta %>%
  group_by(condition) %>%
  summarise(total_cells = n(), .groups = "drop")

# ==========================================================
# Calculate percentages
# ==========================================================

pct_data <- count_data %>%
  left_join(total_cells, by = "condition") %>%
  mutate(
    percentage = 100 * cell_count / total_cells
  )

# ==========================================================
# Create CSV table
# ==========================================================

csv_output <- pct_data %>%
  select(cell_type, condition, percentage) %>%
  pivot_wider(
    names_from = condition,
    values_from = percentage,
    values_fill = 0
  ) %>%
  rename(
    Healthy_Percent = Normal,
    Disease_Percent = Disease
  ) %>%
  mutate(
    Healthy_Percent = round(Healthy_Percent, 2),
    Disease_Percent = round(Disease_Percent, 2)
  )

write.csv(
  csv_output,
  file.path(result_dir,
            "CellType_Percentage_of_Total_Healthy_and_Disease_Cells.csv"),
  row.names = FALSE
)

# ==========================================================
# Plot 1: Side-by-side percentage comparison
# ==========================================================

plot_data <- pct_data

p1 <- ggplot(
  plot_data,
  aes(
    x = cell_type,
    y = percentage,
    fill = condition
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8)
  ) +
  labs(
    title = "Cell Type Composition by Condition",
    x = "Cell Type",
    y = "Percentage of Total Cells",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# ==========================================================
# Plot 2: Stacked composition
# ==========================================================

p2 <- ggplot(
  meta,
  aes(
    x = condition,
    fill = cell_type
  )
) +
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    title = "Cell Type Composition Within Each Condition",
    x = "Condition",
    y = "Percentage",
    fill = "Cell Type"
  ) +
  theme_minimal()

# ==========================================================
# Save plots
# ==========================================================

pdf(
  file.path(result_dir,
            "CellType_Composition_Plots.pdf"),
  width = 12,
  height = 7
)

print(p1)
print(p2)

dev.off()

cat(format(Sys.time(), "[%H:%M:%S]"), "Done!\n")