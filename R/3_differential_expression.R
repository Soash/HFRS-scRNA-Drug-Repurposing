library(Seurat)
library(dplyr)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/HTNV")

result_dir <- "R/result/3_de"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 3: DIFFERENTIAL EXPRESSION (DISEASE VS NORMAL)
### ==================================================================

cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING PHASE 3: DIFFERENTIAL EXPRESSION\n")

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading Integrated Object...\n")
obj <- qs_read("R/result/2_integration/2_integration.qs2")

# ASGARD requires DE on the raw RNA counts, not the integrated assay
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj) # Ensure layers are joined for FindMarkers

# Create a combined metadata column for easy subsetting (e.g., "0_Disease")
obj$cluster_condition <- paste0(obj$seurat_clusters, "_", obj$condition)
Idents(obj) <- "cluster_condition"

# Get a list of all unique clusters
clusters <- sort(unique(as.character(obj$seurat_clusters)))

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step 2: Calculating DEGs per cluster (Disease vs Normal)...\n")

for (cluster_id in clusters) {
  
  ident_disease <- paste0(cluster_id, "_Disease")
  ident_normal <- paste0(cluster_id, "_Normal")
  
  # Check if both conditions exist in this cluster (some rare cell types might only exist in one)
  if (ident_disease %in% unique(Idents(obj)) && ident_normal %in% unique(Idents(obj))) {
    
    cat(format(Sys.time(), "[%H:%M:%S]"), "Processing Cluster", cluster_id, "...\n")
    
    # Run standard Seurat Wilcoxon test (supported by ASGARD)
    # Using a loose logfc.threshold to ensure enough background genes for the ASGARD K-S test
    markers <- FindMarkers(obj, 
                           ident.1 = ident_disease, 
                           ident.2 = ident_normal, 
                           logfc.threshold = 0.1, 
                           min.pct = 0.1,
                           verbose = FALSE)
    
    # Format for ASGARD: It typically expects the gene names as a column, not just rownames
    markers$Gene <- rownames(markers)
    
    # Save the DE table
    output_file <- file.path(result_dir, paste0("Cluster_", cluster_id, "_DEGs.csv"))
    write.csv(markers, file = output_file, row.names = FALSE)
    
  } else {
    cat(format(Sys.time(), "[%H:%M:%S]"), "Skipping Cluster", cluster_id, "- missing cells in one condition.\n")
  }
}

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "PHASE 3 COMPLETE! DE files saved to", result_dir, "\n")

sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()


