library(tidyverse)
library(Seurat) # technically do not need this but 
library(ggrepel)


#' scRNA-seq Volcano Plot
#'
#' This function creates a volcano plot for differential expression results.
#' Input can be the resulting dataframe from `FindAllMarkers` or a saved csv file which is read in using `read_csv`.
#' Make sure there is column with gene names labeled as "gene".
#'
#' Function depends on ggplot and ggrepel.
#'
#' @param deg_df A dataframe containing results from differential expression
#' @param lfc_cutoff A number. Cutoff for highlighting points based on log2FoldChange. Default=1.5
#' @param adj_pval_cutoff A number. Cutoff for highlighting points based on -log10(adjusted p-values). Default=2
#' @param title A string of text. Here you can specify a plot title. Default=NULL
#' @param st A string of text. Here you can specify a plot subtitle. Default=NULL
#' @param cap A string of text. Here you can specify a plot caption. Default=NULL
#'
#' @return Volcano plot.
#' @export
#'
#' @examples
#' scrnaseq_vp(deg_dataframe, title='Hour 1', subtitle='Celltype A vs Celltype B', caption='(Left) Celltype B || Celltype A (Right)')
#' 
#' 

scrnaseq_vp <- function(deg_df, lfc_cutoff=1.5, adj_pval_cutoff = 2, title=NULL, st=NULL, cap=NULL){
  
  vp <- ggplot(deg_df, aes(avg_log2FC, -log10(p_val_adj))) +
    geom_point(color = 'gray') +
    geom_point(data = subset(deg_df, avg_log2FC > lfc_cutoff & -log10(p_val_adj) > adj_pval_cutoff), 
               aes(avg_log2FC, -log10(p_val_adj)), color = 'black') + 
    ggrepel::geom_text_repel(data = subset(deg_df, avg_log2FC > lfc_cutoff & -log10(p_val_adj) > adj_pval_cutoff), 
                             aes(avg_log2FC, -log10(p_val_adj), label = gene)) + 
    geom_point(data = subset(deg_df, avg_log2FC < -lfc_cutoff & -log10(p_val_adj) > adj_pval_cutoff), 
               aes(avg_log2FC, -log10(p_val_adj)), color = 'black') + 
    ggrepel::geom_text_repel(data = subset(deg_df, avg_log2FC < -lfc_cutoff & -log10(p_val_adj) > adj_pval_cutoff), 
                             aes(avg_log2FC, -log10(p_val_adj), label = gene)) +
    geom_hline(yintercept = sig_p, color = 'red', linetype = 'dashed') +
    labs(title = title, 
         subtitle = st,
         caption = cap) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text = element_text(color = 'black'))
  return(vp)
}