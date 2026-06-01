# ==================================================================
# 1. Environment Setup
# ==================================================================
rm(list = ls())
gc(); gc()

library(Seurat)
library(Asgard) # Corrected capitalization
library(dplyr)
library(qs2)

setwd("~/HTNV")

result_dir <- "R/result/4_Asgard"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### ASGARD PIPELINE: HFRS DRUG REPURPOSING
### ==================================================================
cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING ASGARD PIPELINE FOR HFRS\n")

# --- Step 1: Load Data & Prepare Assay ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 1: Loading integrated Seurat object...\n")
combined <- qs_read("R/result/2_integration/2_integration.qs2")

# ASGARD requires standard RNA normalization for FindMarkers to work correctly
DefaultAssay(combined) <- "RNA"

cat(format(Sys.time(), "[%H:%M:%S]"), "Joining layers for Differential Expression...\n")
combined <- JoinLayers(combined)
combined <- NormalizeData(combined, verbose = FALSE)

# --- Step 2: Calculate Cell-Type Specific DEGs (Disease vs Normal) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 2: Calculating DEGs per Seurat Cluster...\n")

Gene.list <- list()
C_names <- NULL
clusters <- sort(unique(as.character(combined$seurat_clusters)))

for(i in clusters){
  Idents(combined) <- "seurat_clusters"
  c_cells <- subset(combined, idents = i)
  
  counts <- table(c_cells$condition)
  
  if("Disease" %in% names(counts) && "Normal" %in% names(counts) && counts["Disease"] > 10 && counts["Normal"] > 10){
    
    cat("Processing Cluster:", i, "... ")
    Idents(c_cells) <- "condition"
    
    # Standard Calculation: Disease vs Normal within the same cluster
    C_data <- FindMarkers(c_cells, ident.1 = "Disease", ident.2 = "Normal", verbose = FALSE)
    
    # Format for ASGARD requirements
    Gene.list[[i]] <- data.frame(
      row.names = row.names(C_data),
      score = C_data$avg_log2FC,
      adj.P.Val = C_data$p_val_adj,
      P.Value = C_data$p_val
    )
    C_names <- c(C_names, i)
    cat("Done.\n")
    
  } else if (i == "17") {
    
    # ===============================================================
    # THE FIX: CLUSTER 17 (NRBCs) SPECIAL EXCEPTION
    # Because healthy patients have 0 NRBCs, we compare Disease NRBCs 
    # against ALL Normal PBMCs to capture the abnormal signature.
    # ===============================================================
    cat("Processing Cluster:", i, "(NRBC Special Case)... ")
    
    # Set the identity on the global object to condition
    Idents(combined) <- "condition"
    
    # Isolate barcodes: Disease cells that belong to Cluster 17
    nrbc_cells <- intersect(WhichCells(combined, idents = "Disease"), 
                            WhichCells(combined, expression = seurat_clusters == i))
    
    # Isolate barcodes: All Normal baseline cells
    normal_cells <- WhichCells(combined, idents = "Normal")
    
    # Run DE: Disease NRBCs vs All Normal Cells
    C_data <- FindMarkers(combined, ident.1 = nrbc_cells, ident.2 = normal_cells, verbose = FALSE)
    
    # Format for ASGARD requirements
    Gene.list[[i]] <- data.frame(
      row.names = row.names(C_data),
      score = C_data$avg_log2FC,
      adj.P.Val = C_data$p_val_adj,
      P.Value = C_data$p_val
    )
    C_names <- c(C_names, i)
    cat("Done.\n")
    
  } else {
    cat("Skipping Cluster", i, "- Insufficient cells in one or both groups.\n")
  }
}
names(Gene.list) <- C_names

# --- Step 3: Drug Repurposing (Mono-Drug) ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 3: Loading Drug Reference & Running Repurposing...\n")

# Use the exact LINCS haematopoietic reference files found in your directory
my_gene_info <- read.table("database/L1000/DrugReference/haematopoietic-and-lymphoid-tissue_gene_info.txt", sep="\t", header = TRUE, quote = "")
my_drug_info <- read.table("database/L1000/DrugReference/haematopoietic-and-lymphoid-tissue_drug_info.txt", sep="\t", header = TRUE, quote = "")

drug.ref.profiles <- GetDrugRef(
  drug.response.path = 'database/L1000/DrugReference/haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
  probe.to.genes = my_gene_info, 
  drug.info = my_drug_info
)

Drug.ident.res <- GetDrug(
  gene.data = Gene.list, 
  drug.ref.profiles = drug.ref.profiles, 
  repurposing.unit = "drug", 
  connectivity = "negative", 
  drug.type = "FDA"
)

qs_save(Drug.ident.res, file = file.path(result_dir, "Drug.ident.res.qs2"), nthreads = 12)















# --- Step 4: Final Scoring & Fixed Export ---
cat(format(Sys.time(), "[%H:%M:%S]"), "Step 4: Calculating Therapeutic Scores & Saving CSV...\n")

cell_metadata <- combined@meta.data
# Map clusters and samples to the metadata structure ASGARD expects
cell_metadata$cluster <- as.character(combined$seurat_clusters)
cell_metadata$sample <- combined$orig.ident 

# Define which samples are the diseased cases to anchor the therapeutic scoring
Case_samples <- unique(cell_metadata$sample[cell_metadata$condition == "Disease"])

# ===============================================================
# THE FIX: DYNAMIC TISSUE STRING DETECTION
# Auto-detect the exact LINCS spelling (spaces vs underscores)
# ===============================================================
cinfo_path <- "database/L1000/GSE92742_Broad_LINCS_cell_info.txt"
cinfo <- read.delim(cinfo_path, stringsAsFactors = FALSE)

# Look for the native hematopoietic/lymphoid string used in the file
valid_tissues <- unique(cinfo$primary_site)
lincs_tissue_str <- valid_tissues[grep("haematopoietic|hematopoietic|lymphoid|blood", valid_tissues, ignore.case = TRUE)][1]

cat(format(Sys.time(), "[%H:%M:%S]"), "Auto-detected LINCS primary site:", lincs_tissue_str, "\n")

# Run DrugScore with ONLY the accepted arguments
Drug.score <- DrugScore(
  cell_metadata, 
  cluster_degs = Gene.list, 
  cluster_drugs = Drug.ident.res, 
  tissue = lincs_tissue_str, # Use the dynamically detected perfect string
  case = Case_samples, 
  gse92742_gctx_path = "database/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
  gse70138_gctx_path = "database/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
)

# Convert row names to a permanent column and filter for significant therapeutic candidates
Final.drugs <- Drug.score %>%
  mutate(Drug_Name = rownames(.)) %>%
  filter(FDR < 0.05) %>%
  arrange(desc(Drug.therapeutic.score)) %>%
  dplyr::select(Drug_Name, Drug.therapeutic.score, P.value, FDR)

# Save the clean CSV
write.csv(Final.drugs, file.path(result_dir, "ASGARD_HFRS_Candidates.csv"), row.names = FALSE)

cat(format(Sys.time(), "[%H:%M:%S]"), "ASGARD PIPELINE COMPLETE!\n")

# --- Teardown ---
sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()







