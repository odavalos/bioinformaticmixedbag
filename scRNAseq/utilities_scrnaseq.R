# Useful scRNA-seq Functions

#' Number of unique clusters generated per resolution
#' 
#' @param seurat_data A Seurat object
#' @param min_res A logical indicating whether to identify resolution with minimum number of clusters
#' 
#' @return Returns either the size of unique clusters generated or specific resolution with smallest unique cluster size
#' 
#' @example 
#' getuniqueclustersizes(seurat_data, min_res=FALSE)

getuniqueclustersizes <- function(seurat_data, min_res = FALSE) {
  if (min_res == TRUE) {
    # find number of unique clusters for each cluster resolution
    res <-
      grep("res", colnames(seurat_data@meta.data), value = TRUE)
    cluster_nums <- res %>% purrr::map_chr( ~ length(unique(seurat_data@meta.data[, .x])))
    
    # make a dataframe for unique clusters
    df <- data.frame(res = res, clusters = cluster_nums)
    
    # get index for smallest resolution
    idx <- which.min(df$clusters)
    
    # pull the unique resolution
    min_r <- df[idx, 'res']
    
    return(min_r)
    
  } else if (min_res = FALSE)
    
    # grep all columns that contain 'res' in metadata and determin the size of clusters generated
    resolutions <-
      grep("res", colnames(seurat_data@meta.data), value = TRUE) %>%
      purrr::map_chr( ~ paste(.x, "--> clusters generated:", length(unique(
        seurat_data@meta.data[, .x]
      ))))
  return(resolutions)
} 



#' Get percentage of gene expression greater than zero in dataset
#'
#' @param seurat_obj A Seurat object
#' @param gene A gene of interest
#' 
#' @return A plot containing a UMAP from the seurat object and either density or expression plots
#' 
#' 
#' @example 
#' getgenetotals(seurat_data, gene)
#' 


getgenetotals <- function(seurat_data, gene) {
  # evaluate total population of gene+ cells in subclusters
  
  # Get total number of cells in the seurat dataset
  total_cells <- seurat_data %>%
    ncol()
  
  # Calculate total number cells in which contain expression greater than zero in subclusters
  hits <- 0 # hits serves as our counter start
  gene_cells <-
    ifelse(seurat_data[gene]@assays$RNA@counts > 0, hits + 1, hits + 0) %>% sum() # counter
  
  # Print message on the console with percent of gene found in subclusters
  message(
    'Percent of cells with \'',
    gene,
    '\' expression > 0 in subclusters : ',
    round(gene_cells / total_cells * 100, 2),
    "%"
  )
  
}



# Plotting Functions ------------------------------------------------------


#' Plot gene expression for enhanced visualization
#'
#' @param seurat_obj A Seurat object
#' @param genelist A list of genes in vector form
#' @param density A logical indicating whether to use density plots
#' 
#' @return A plot containing a UMAP from the seurat object and either density or expression plots
#' 
#' 
#' @examples 
#' genes <- c('gene1', 'gene2', 'gene3')
#' plot_GexViz(seurat_obj, genes, density = TRUE)
#' 
#' gene2
#' plot_GexViz(seurat_obj, 'gene2', density = FALSE)

plot_GexViz <- function(seurat_obj, genelist, density = FALSE) {
  # visualize gene expression for annotations using umap of clusters and density or expression plots
  
  require(ggthemes)
  require(Nebulosa)
  
  
  if (density == TRUE) {
    # plot clusters and gene density
    
    p_umap <- DimPlot(
      seurat_obj,
      reduction = "umap",
      label = TRUE,
      repel = TRUE,
      label.size = 8
    ) +
      scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
      theme(plot.title = element_text(hjust = 0.5)) +
      NoLegend()
    
    # Density plots
    p_density <- plot_density(seurat_obj, genelist, pal = 'magma') 
    
    # using patchwork `|` output both plots in one
    return(p_umap | p_density)
    
  } else if (density == FALSE) {
    # plot clusters and gene expression
    
    p_umap <- DimPlot(
      seurat_obj,
      reduction = "umap",
      label = TRUE,
      repel = TRUE,
      label.size = 8
    ) +
      scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
      theme(plot.title = element_text(hjust = 0.5)) +
      NoLegend()
    
    # Expression plot 
    p_gex <- FeaturePlot(seurat_obj, features = genelist) 
    
    # using patchwork `|` output both plots in one
    return(p_umap | p_gex)
    
  }
}




#' Plot cluster proportions
#'
#' @param seurat_obj A Seurat object
#' @param clust_res A clustering result object
#' @param annotated A logical indicating whether seurat object is annotated
#' @param percents A logical indicating whether to plot percentages on barplots
#' @param integrated A logical indicating whether seurat object is integrated
#' 
#' @return A plot or plots of cluster proportions in the input seurat object
#' @export
#'
#'
#' @examples
#' plotclusterprops(seurat_obj, clust_res)
#' plotclusterprops(seurat_obj, clust_res, annotated=TRUE)
#' plotclusterprops(seurat_obj, clust_res, percents=TRUE)
#' plotclusterprops(seurat_obj, clust_res, integrated=TRUE)
#' plotclusterprops(seurat_obj, clust_res, annotated=TRUE, percents=TRUE)
#' plotclusterprops(seurat_obj, clust_res, annotated=TRUE, integrated=TRUE)
#' plotclusterprops(seurat_obj, clust_res, annotated=TRUE, percents=TRUE, integrated=TRUE)


plotclusterprops <-  function(seurat_obj, clust_res, annotated=FALSE, percents=FALSE, integrated=FALSE) { 
  
  ## take an integrated Seurat object, plot distributions over orig.ident
  require(RColorBrewer)
  
  # Integrated Plots
  if(integrated==TRUE){
    
    # make plot for annotated clusters  
    if(annotated==TRUE){
      
      # get cluster proportions
      merged_mat <- seurat_obj@meta.data %>%
        group_by(celltypes, orig.ident) %>%
        select(orig.ident, celltypes) %>%
        summarize(value = n())
      
      colnames(merged_mat)[2] <- "dataset"
      
      # order by dataset
      merged_mat2 <- merged_mat %>%
        extract(dataset, c("dataset_num"), regex = '([0-9]+)', remove = FALSE, convert = TRUE) %>%
        arrange(dataset_num, desc(value)) %>%
        select(-dataset_num)
      
      # pull out ordered names
      lvs <- unique(merged_mat2$dataset)
      
      # convert dataset into factor and relevel it
      merged_mat2  <- merged_mat2 %>%
        mutate(dataset = as.factor(dataset),
               dataset = fct_relevel(dataset, lvs))
      
      # calculate percentages
      merged_mat2 <- merged_mat2 %>%
        group_by(dataset) %>%
        mutate(percent=value/sum(value))
      
      # arrange clusters by size
      merged_mat3 <- merged_mat %>% 
        ungroup() %>% 
        group_by(celltypes) %>% 
        summarize(value = sum(value)) %>% 
        arrange(desc(value))
      
      # stacked barplot with percentages
      if(percents==TRUE){
        
        p3 <- ggplot(merged_mat2,aes(x=reorder(celltypes,value),y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity", 
                   color = 'black', alpha = 0.9) + 
          theme_bw() + 
          coord_flip() + 
          scale_fill_brewer(palette = "Dark2", direction = 1) +
          labs(x="Cluster", y="Cell fraction per dataset",
               fill = "") +
          theme(legend.position="top",
                axis.text = element_text(color = 'black'))
        p3 <- p3 + geom_text(aes(label=paste0(sprintf("%1.1f", percent*100),"%")),
                             position=position_fill(vjust=0.5), colour="white", size = 3)
        
        # stacked barplot (no percentages)  
      } else if(percents==FALSE){
        
        p3 <- ggplot(merged_mat2,aes(x=reorder(celltypes,value),y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity", 
                   color = 'black', alpha = 0.9) + 
          theme_bw() + 
          coord_flip() + 
          scale_fill_brewer(palette = "Dark2", direction = 1) +
          labs(x="Cluster", y="Cell fraction per dataset",
               fill = "") +
          theme(legend.position="top",
                axis.text = element_text(color = 'black'))
        
      }
      
      # cell per cluster dotplot
      p4 <- ggplot(merged_mat3, aes(y= reorder(celltypes, value),x = value, fill = value, size = value)) + 
        geom_point(shape=21) + 
        theme_bw() + 
        scale_x_log10() + 
        labs(x="Cells per cluster (log10)", y="",
             fill = "Cells", size="Cells") +
        scale_fill_viridis_c(option = 'viridis') +
        theme(axis.text = element_text(color = 'black'))
      
      p3 + p4 + patchwork::plot_layout(widths = c(3,1))
      
      # make plot for unannotated clusters  
    } else if(annotated==FALSE){
      
      # get cluster proportions
      cts_tbl <- table(seurat_obj@meta.data[,chosen_res], 
                       seurat_obj@meta.data$orig.ident)
      cts_mat <- as.data.frame.matrix(cts_tbl)
      cts_mat$cluster <- rownames(cts_mat)
      merged_mat <- cts_mat %>%
        pivot_longer(-cluster, names_to = 'variable', values_to = 'value')
      merged_mat$cluster <- as.factor(merged_mat$cluster)
      
      # get cluster sizes
      cluster_size <- merged_mat %>%
        group_by(cluster) %>%
        summarise(value = sum(value))
      # order by size
      sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
      cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
      merged_mat$cluster <- factor(merged_mat$cluster,levels = sorted_labels)
      
      colnames(merged_mat)[2] <- "dataset"
      
      # order by dataset
      merged_mat <- merged_mat %>%
        extract(dataset, c("dataset_num"), regex = '([0-9]+)', remove = FALSE, convert = TRUE) %>%
        arrange(dataset_num, desc=FALSE) %>%
        select(-dataset_num)
      
      # pull out ordered names
      lvs <- unique(merged_mat$dataset)
      
      # convert dataset into factor and relevel it
      merged_mat  <- merged_mat %>%
        mutate(dataset = as.factor(dataset),
               dataset = fct_relevel(dataset, lvs))
      
      # calculate percentages
      merged_mat <- merged_mat %>%
        group_by(dataset) %>%
        mutate(percent=value/sum(value))
      
      # stacked barplot
      if(percents==TRUE){
        
        p1 <- ggplot(merged_mat,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity", 
                   color = 'black', alpha = 0.9) + 
          theme_bw() + 
          coord_flip() + 
          scale_fill_brewer(palette = "Dark2", direction = 1) +
          labs(x="Cluster number", y="Cell fraction per dataset",
               fill = "") +
          theme(legend.position="top",
                axis.text = element_text(color = 'black'))
        p1 <- p1 + geom_text(aes(label=paste0(sprintf("%1.1f", percent*100),"%")),
                             position=position_fill(vjust=0.5), colour="white", size = 3)
        
        # stacked barplot (no percentages)
      } else if(percents==FALSE){
        
        p1 <- ggplot(merged_mat,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity", 
                   color = 'black', alpha = 0.9) + 
          theme_bw() + 
          coord_flip() + 
          scale_fill_brewer(palette = "Dark2", direction = 1) +
          labs(x="Cluster number", y="Cell fraction per dataset",
               fill = "") +
          theme(legend.position="top",
                axis.text = element_text(color = 'black'))
        
      }
      
      
      # cell per clust dotplot
      p2 <- ggplot(cluster_size, aes(y= cluster,x = value, fill = value, size = value)) + 
        geom_point(shape=21) + 
        theme_bw() + 
        scale_x_log10() + 
        labs(x="Cells per cluster (log10)", y="",
             fill = "Cells", size="Cells") +
        
        scale_fill_viridis_c(option = 'viridis') +
        theme(axis.text = element_text(color = 'black'))
      
      p1 + p2 + patchwork::plot_layout(widths = c(3,1))
      
    }
    
    # plots for non integrated datasets
  } else if(integrated==FALSE){
    
    # calculate celltype proportions
    celltype_props <- data.filt@meta.data %>%
      group_by(celltypes) %>%
      summarize(value = n()) %>%
      arrange(value) %>%
      mutate(percent=value/sum(value))
    
    # sort by celltype and size
    sorted_labels <- celltype_props %>% pull(celltypes)
    celltype_props$celltypes <- factor(celltype_props$celltypes,levels = sorted_labels)
    
    # make dotplot
    dplot <- ggplot(celltype_props, aes(y= celltypes,x = value, fill = value, size = value)) + 
      geom_point(shape=21) + 
      theme_bw(base_size = 14) + 
      scale_x_log10() + 
      labs(x="Cells per cluster (log10)", y="",
           fill = "Cells", size="Cells") +
      
      scale_fill_viridis_c(option = 'viridis') +
      theme(axis.text = element_text(color = 'black'))
    
    dplot
    
  }
  
}
