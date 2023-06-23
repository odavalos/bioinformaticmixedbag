#! /usr/bin/env Rscript
# 
# This script is for running TradeSeq on Seurat and Slingshot objects
# to generate the GAM models.
# The script will run evaluateK and fitGAM. Depending on the arguments.
#
# Usage: Rscript run_tradeSeq.r -s seurat_obj -l slingshot_obj -n samplename -e evaluateK -g fitGAM -k knots -o out_path
#
# Arguments:
#   -s, --seurat_obj: path to seurat object
#   -l, --slingshot_obj: path to slingshot object
#   -n, --samplename: name of sample
#   -e, --evaluateK: whether to run evaluateK
#   -g, --fitGAM: whether to run fitGAM
#   -k, --knots: number of knots for GAM
#   -o, --out_path: path to output directory
#
# Output:
#   - samplename_tradeseq_evaluateK.pdf: plot of evaluateK results
#   - samplename_tradeseq_GAM_sce.rds: sce object with GAM fit
#   TODO: - samplename_tradeseq_console_output.txt: console output 
#
# Dependencies:
#   - R
#   - tidyverse (not really needed)
#   - optparse
#   - Seurat
#   - slingshot
#   - tradeSeq


# Load libraries
library(tidyverse)
library(optparse)
library(Seurat)
library(slingshot)
library(tradeSeq)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser, 
                     c("--seurat_obj","-s"), 
                     type="character", 
                     default="data.rds")
parser <- add_option(parser,
                        c("--slingshot_obj","-l"),
                        type="character",
                        default="data.rds")
parser <- add_option(parser,
                        c("--samplename","-n"),
                        type="character",
                        default="sample")
parser <- add_option(parser,
                        c("--evaluateK","-e"),
                        type="logical",
                        default=TRUE)
parser <- add_option(parser,
                        c("--fitGAM","-g"),
                        type="logical",
                        default=TRUE)
parser <- add_option(parser,
                        c("--knots","-k"),
                        type="integer",
                        default=6)
parser <- add_option(parser,
                     c("--out_path","-o"),
                     type="character",
                     default="out.rds")

# Parse the command-line arguments.
out    <- parse_args(parser)
seurat_obj <- out$seurat_obj
slingshot_obj <- out$slingshot_obj
samplename <- out$samplename
evaluateK <- out$evaluateK
GAM <- out$fitGAM
knots <- out$knots
out_p <- out$out_path


# print out the arguments
cat("\n")
cat("seurat_obj: ", seurat_obj, "\n")
cat("slingshot_obj: ", slingshot_obj, "\n")
cat("samplename: ", samplename, "\n")
cat("out: ", out_p, "\n")
cat("\n")


# load seurat object
seurat <- readRDS(seurat_obj)

# load slingshot object
sds <- readRDS(slingshot_obj)

# get the log transformed counts matrix
seurat_cts <- GetAssayData(seurat, slot = "data")[VariableFeatures(seurat),]

# set the seed for reproducibility
set.seed(2023)

# set up parallelization
BPPARAM <- BiocParallel::register(BiocParallel::MulticoreParam())

# run tradeSeq's evaluateK

if(evaluateK == TRUE){


  pdf(file = paste0(out_p, samplename, "_tradeseq_evaluateK.pdf"), height = 8, width = 12)
  icMat <- evaluateK(
    counts = seurat_cts,
    sds = sds,
    k = 3:20,
    nGenes = 200,
    verbose = T,
    parallel = T,
    BPPARAM = BPPARAM$MulticoreParam
  )
    dev.off()

} 


# run tradeSeq's fitGAM

if(GAM == TRUE){
    sceGAM <- fitGAM(
        counts = seurat_cts,
        sds = sds,
        nknots = knots,
        verbose = T,
        parallel = T,
        sce = TRUE,
        BPPARAM = BPPARAM$MulticoreParam
    )

    # save gam sce object
    saveRDS(sceGAM, file = paste0(out_p, samplename, "_tradeseq_GAM_sce.rds"))

    }

print("Done!")


# # save all console output
# sink(file = paste0(out, "_", samplename, "_tradeseq_console_output.txt"), append = TRUE)
