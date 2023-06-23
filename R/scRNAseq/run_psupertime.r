#! /usr/bin/env Rscript
#
# This script is for running psupertime on Seurat objects
# to generate a supervised pseudotime for a celltype.
#
# Usage: Rscript run_psupertime.r -s seurat_obj -n samplename -c celltype -o out_path
#
# Arguments:
#   -s, --seurat_obj: path to seurat object
#   -n, --samplename: name of sample
#   -c, --celltype: name of celltype pseudotime is being generated for
#   -o, --out_path: path to output directory
#
# Output:
#   - samplename_celltype_psupertime.rds: psupertime object
# 
# Dependencies:
#   - R
#   - tidyverse
#   - optparse
#   - Seurat
#   - psupertime (https://github.com/wmacnair/psupertime)
#



library(tidyverse)
library(optparse)
library(Seurat)
library(psupertime)


# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser, 
                     c("--seurat_obj","-s"), 
                     type="character", 
                     default="data.rds")
parser <- add_option(parser,
                        c("--samplename","-n"),
                        type="character",
                        default="sample")
parser <- add_option(parser,
                        c("--celltype","-c"),
                        type="character",
                        default="celltype")
parser <- add_option(parser,
                     c("--out_path","-o"),
                     type="character",
                     default="out.rds")

# Parse the command-line arguments.
out    <- parse_args(parser)
seurat_obj <- out$seurat_obj
samplename <- out$samplename
celltype <- out$celltype
out_p <- out$out_path


# print out the arguments
cat("\n")
cat("seurat_obj: ", seurat_obj, "\n")
cat("samplename: ", samplename, "\n")
cat("celltype: ", celltype, "\n")
cat("out_path: ", out_p, "\n")
cat("\n")


# load seurat object
seurat <- readRDS(seurat_obj)

# extract variable genes
vargenes <- VariableFeatures(seurat)

# convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat, assay = 'RNA')

# subset to variable genes
sce <- sce[vargenes,]

# run psupertime
y <- sce$Timepoint # already a factor
psuper_obj <- psupertime(sce, y, sel_genes="all", seed = 2023)

# save psupertime  object
saveRDS(psuper_obj, file = paste0(out_p, samplename, "_", celltype, "_psupertime.rds"))

print("Done!")
