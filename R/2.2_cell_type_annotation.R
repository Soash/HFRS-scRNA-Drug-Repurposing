library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)

# ==================================================================
# Environment Setup
# ==================================================================
rm(list = ls())
gc(); gc()

setwd("~/HTNV")

# Create new directory for this phase
result_dir <- "R/result/2.2_annotation"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 2.2: CLUSTER IDENTIFICATION & CELL TYPE ANNOTATION
### ==================================================================

cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 2.2: CELL TYPE ANNOTATION\n")

# --- Step 1: Load Integrated Data ---
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading Integrated Object...\n")
combined <- qs_read("R/result/2_integration/2_integration.qs2")

# Ensure we are working with standard RNA counts for marker identification
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)

# --- Step 2: Find All Cluster Markers ---
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 2: Calculating global markers for all clusters...\n")
# Using the criteria outlined in the HFRS paper: lnFC > 0.25, min.pct > 0.1, wilcox test
cluster_markers <- FindAllMarkers(combined, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25,
                                  test.use = "wilcox",
                                  verbose = TRUE)

# Save the full marker list for your records
write.csv(cluster_markers, file.path(result_dir, "All_Cluster_Markers.csv"), row.names = FALSE)

# Extract and save the top 10 markers per cluster for easy manual review
top10_markers <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10_markers, file.path(result_dir, "Top10_Cluster_Markers.csv"), row.names = FALSE)

# --- Step 3: Map Numerical Clusters to Biological Cell Types ---
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 3: Applying Cell Type Annotations...\n")

cluster_to_celltype <- c(
  "0" = "CD8+ T cells",
  "1" = "Naive T cells",
  "2" = "Cycling cells",
  "3" = "Cycling cells",
  "4" = "Macrophages",
  "5" = "Classical Monocytes",
  "6" = "NK cells",
  "7" = "Plasma cells",
  "8" = "Monocytes",
  "9" = "B cells",
  "10" = "Naive T cells",
  "11" = "T cells",
  "12" = "Gamma-delta T cells",
  "13" = "Naive B cells",
  "14" = "NK cells",
  "15" = "Dendritic cells",
  "16" = "Dendritic cells",
  "17" = "NRBCs",
  "18" = "Megakaryocytes",
  "19" = "Basophils"
)

# Rename the identities in the Seurat object
names(cluster_to_celltype) <- levels(combined)
combined <- RenameIdents(combined, cluster_to_celltype)

# Save this explicitly to the metadata so it persists in downstream scripts
combined$cell_type <- Idents(combined)

# --- Step 4: Generate Annotated Plots ---
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 4: Generating Annotated UMAPs...\n")

pdf(file.path(plot_dir, "Annotated_UMAP.pdf"), width = 10, height = 8)
print(DimPlot(combined, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 5, repel = TRUE) + 
        ggtitle("UMAP by Annotated Cell Type"))
dev.off()

# Plot canonical markers to verify your annotations are correct
canonical_markers <- c("CD3E", "CD14", "MS4A1", "MZB1", "NKG7", "PPBP", "HBA1", "HBB")
pdf(file.path(plot_dir, "Verification_FeaturePlots.pdf"), width = 16, height = 8)
print(FeaturePlot(combined, features = canonical_markers, reduction = "umap", ncol = 4, pt.size = 0.5))
dev.off()

# --- Step 5: Save Annotated Object ---
cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 5: Saving Annotated Object...\n")
qs_save(combined, file = file.path(result_dir, "2.2_annotation.qs2"), nthreads = 12)

cat(format(Sys.time(), "[%H:%M:%S]"), "PHASE 2.2 COMPLETE!\n")

# --- Teardown ---
sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()