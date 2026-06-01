library(Seurat)
library(ggplot2)
library(patchwork) # Handles side-by-side layouts
library(qs2)

# Set working directory
setwd("~/HTNV")

# Define plot output path
output_png <- "R/result/2.2_annotation/plots/Phase2_Integration_UMAPs.png"

# 1. Load the annotated Seurat object
combined <- qs_read("R/result/2.2_annotation/2.2_annotation.qs2")

# 2. Generate the three individual plots
p1 <- DimPlot(combined, group.by = "orig.ident") + 
  ggtitle("UMAP by Patient") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- DimPlot(combined, group.by = "condition") + 
  scale_color_manual(values = c("Normal" = "#4DAF4A", "Disease" = "#E41A1C")) +
  ggtitle("UMAP by Condition") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p3 <- DimPlot(combined, group.by = "cell_type", label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle("UMAP by Cell Type") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

# 3. Open the PNG device with a wide layout aspect ratio (3:1)
png(output_png, width = 18, height = 6, units = "in", res = 300)

# 4. Use patchwork syntax '|' to place them side by side
print(p1 | p2 | p3)

# 5. Close and save the file
dev.off()

cat("Successfully saved side-by-side UMAP matrix to:", output_png, "\n")