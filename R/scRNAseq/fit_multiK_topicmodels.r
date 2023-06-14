#! /usr/bin/env Rscript
# 
# Adapted from KLRhodes's https://github.com/KLRhodes script `prefit_poisson_nmf.R` at 
# https://github.com/KLRhodes/Embryoid_Body_Pilot_Workflowr/blob/master/code/prefit_poisson_nmf.R
# 
# This script fits topic models using fastTopics and iterates
# over a range of k values and saves the outputs. Need to run on 
# a cluster since I keep getting this error: 
# `Quick-TRANSfer stage steps exceeded maximum (= 804250)`. Once 
# the models are fit, the models can be loaded and analyzed in
# on a local machine.
#


# Load a few packages.
library(Matrix)
library(optparse)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser, 
                     c("--input","-i"), 
                     type="character", 
                     default="data.rds")
parser <- add_option(parser,
                     c("--out","-o"),
                     type="character",
                     default="out.rds")
# parser <- add_option(parser,c("--k","-k"),type = "integer",default = 3)
parser <- add_option(parser,
                     c("--method","-m"),
                     type = "character",
                     default = "em")
parser <- add_option(parser,
                     c("--num_iter","-n"),
                     type = "integer",
                     default = 1000)
parser <- add_option(parser,
                     "--nc",
                     type = "integer",
                     default = 4)

# Parse the command-line arguments.
out    <- parse_args(parser)
input       <- out$input
method      <- out$method
num_iter    <- out$num_iter
outfile     <- out$out
# k           <- out$k
nc          <- out$nc

# cat the command-line arguments.
cat("--------------------------------------------\n")
cat("Command-line arguments:\n")
cat("Input file: ", input, "\n")
cat("Output file: ", outfile, "\n")
cat("Method: ", method, "\n")
cat("Number of iterations: ", num_iter, "\n")
cat("Number of cores: ", nc, "\n")
cat("--------------------------------------------\n")


# Set the random seed.
set.seed(2023)
cat("Random seed set to 2023.\n")

# Load the data.
cat("Loading data.\n")
counts_mat <- readRDS(input)
cat("Loaded", 
    nrow(counts_mat),
    "cells and",
    ncol(counts_mat),
    "genes.\n")

# # check if non-finite values are present. if so, remove them.
# if(any(!is.finite(counts_mat@x)){
#     cat("Non-finite values present. Removing them.\n")
#     counts_mat <- counts_mat[!is.finite(counts_mat@x),]
# } else {
#    cat("No non-finite values present.\n")
# }


# Vector of k values to fit
k_values <- c(4,6,8,10,15,20,25,30)

# Empty list to store the results.
all_fits <- list()

# Measure time it take to each topic model. Then add the time to get total time.
runtime_df <- data.frame(k = integer(), 
                         runtime = numeric())


# Fit the topic models.
for(i in k_values){

    cat("\nFitting topic model with k = ", i, "\n")
    runtime <- system.time(
        fit <- fit_topic_model(counts_mat, 
                               k = i, 
                               numiter.main = num_iter,
                               method.main = method,
                               control.main = list(numiter = 4, nc = nc),
                               verbose = "progressbar")
    )
    # add the runtime to the dataframe
    runtime_df <- rbind(runtime_df, 
                        data.frame(k = i, 
                                   runtime = runtime[3]))
    # save the named fit to the list
    all_fits[[paste0("fit_k_",i)]] <- fit

    cat("Fitting topic model with k = ", i, " took ", runtime[3], " seconds.\n")
}

# Add the total runtime to the dataframe.
runtime_df <- rbind(runtime_df, 
                    data.frame(k = "total", 
                               runtime = sum(runtime_df$runtime)))

cat("Total runtime: ", sum(runtime_df$runtime), " seconds.\n")

# Save the runtime dataframe.
cat("Saving runtime dataframe.\n")
write.csv(runtime_df, 
          file = paste0(outfile,"_topicmodels_runtime.csv"))
          

# Save the results.
cat("Saving results.\n")
saveRDS(all_fits, 
        file = paste0(outfile,"_topicmodels.rds"))


# cat the session info.
sessionInfo()

# Capture all standard output and error.
sink(
    file = paste0(outfile,"_topicmodels.log"), 
    append = FALSE,
    split = FALSE,
    type = "output"
)

