library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(qs2)

# ==================================================================

rm(list = ls())
gc(); gc()

setwd("~/HTNV")

result_dir <- "R/result/1_qc"
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

plot_dir <- file.path(result_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

log_file <- file(file.path(result_dir, "log.txt"), open = "wt")

while (sink.number() > 0) {sink()}
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

### ==================================================================
### PHASE 1: SEQUENTIAL QC (STATIC THRESHOLDS)
### ==================================================================

plan("sequential")

# Get Sample Lists
disease_samples <- list.dirs("GSE161354_RAW/disease", full.names = TRUE, recursive = FALSE)
normal_samples <- list.dirs("GSE161354_RAW/normal", full.names = TRUE, recursive = FALSE)
all_samples <- c(disease_samples, normal_samples)

seurat_list <- list()

# ==================================================================
# BEGIN SEQUENTIAL LOOP
# ==================================================================
for (i in seq_along(all_samples)) {
  
  x <- all_samples[i]
  sample_name <- basename(x)
  
  cat(format(Sys.time(), "[%H:%M:%S]"), "STARTING SAMPLE", i, "of", length(all_samples), ":", sample_name, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step A: Loading 10X data...\n")
  obj <- CreateSeuratObject(counts = Read10X(data.dir = x))
  obj$orig.ident <- sample_name
  obj$condition <- ifelse(grepl("disease", x), "Disease", "Normal")
  
  print(obj)
  before_genes <- sum(rowSums(obj[["RNA"]]$counts) > 0)
  before_cells <- ncol(obj)
  cat("Genes (Raw):", before_genes, "| Cells (Raw):", before_cells, "\n")
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step B: Calculating QC Metrics & Plotting Raw...\n")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- NormalizeData(obj, verbose = FALSE)
  
  raw_pdf <- file.path(plot_dir, paste0("Raw_QC_", sample_name, ".pdf"))
  pdf(raw_pdf, width = 12, height = 6)
  print(
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste(sample_name, "- Raw"))
  )
  dev.off()
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step C: Applying Static Filters (nFeature > 200, MT < 20%)...\n")
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step D: Removing MT genes from matrix...\n")
  non_mt_genes <- grep(pattern = "^MT-", x = rownames(obj), invert = TRUE, value = TRUE)
  obj <- subset(obj, features = non_mt_genes)
  
  print(obj)
  final_genes <- sum(rowSums(obj[["RNA"]]$counts) > 0)
  final_cells <- ncol(obj)
  
  
  cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Step E: Finalizing logs & Plotting Cleaned...\n")
  cat(
    "Sample:", sample_name, "| Condition:", obj$condition[1], "\n",
    "[Thresholds Applied]\n",
    "- Max MT%: 20\n",
    "- Min nFeature: 200\n",
    "[Counts Tracking]\n",
    "- Cells:", before_cells, "(Raw) ->", final_cells, "(Final)\n",
    "- Genes:", before_genes, "(Raw) ->", final_genes, "(Active - MT Removed)\n",
    "--------------------------------------------------\n"
  )
  
  clean_pdf <- file.path(plot_dir, paste0("QC_", sample_name, ".pdf"))
  pdf(clean_pdf, width = 12, height = 6)
  print(
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste(sample_name, "- Cleaned"))
  )
  dev.off()
  
  seurat_list[[sample_name]] <- obj
}

cat("\n", format(Sys.time(), "[%H:%M:%S]"), "Saving QC processed objects list...\n")
qs_save(seurat_list, file = file.path(result_dir, "1_QC_List.qs2"), nthreads = 4)

sink(type = "message")
while (sink.number() > 0) {sink(type = "output")}
close(log_file)

rm(list = ls())
gc(); gc()


