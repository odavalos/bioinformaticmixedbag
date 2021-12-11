## scRNA-seq

Interoperability between Seurat and Scanpy.
 - If going from Seurat to scanpy need to convert Seurat object into AnnData object.




| Seurat       | Scanpy           |
|--------------|------------------|
| nCount_RNA   | n_counts         |
| nFeature_RNA | n_genes          |
| Read10X()    | sc.read_10x_h5() |




Generating mutiple cluster resolutions in both Seurat and Scanpy:


#### Seurat

    # depends on seurat and tidyverse (dplyr & purrr)
    some_dataset <- FindNeighbors(object = some_dataset,
                                    dims = n_dims)

    # run the clustering for the multiple resolutions
    some_dataset <- FindClusters(
      object = some_dataset,
      resolution = seq(0.20,2,0.10) # change resolution range here
    )

    # Number of unique clusters generated per resolution
    grep("res", colnames(some_dataset@meta.data), value = TRUE) %>%
      purrr::map_chr( ~ paste(.x, "--> clusters generated:", length(unique(
        some_dataset@meta.data[, .x]
      ))))

    # set the desired resolution as new identity
    Idents(some_dataset) <- "seurat_snn_res.#.#"




#### Scanpy

Multiple louvain cluster resolution function.


    # Depends on scanpy and numpy.
    def multi_louvain(adata, start_res, stop_res, step_res, plot = False):
        """
        Generate multiple louvain cluster resolutions with one function.
        adata: scRNAseq dataset
        start_res: starting value for resolution sequence
        stop_res: stopping value for resolution sequence
        step_res: value for incremental steps in sequence
        plot: generates plot for all resolutions generated assumes umap dim reduction `sc.pl.umap(adata, color=res_list)`
        """
        # create dictionary
        resolutions = {}
        for i in np.arange(start_res, stop_res, step_res): # change resolutions here
            res = round(i,2) # each resolution value
            name = 'louvain_' + str(round(i,2)) # key
            d = {name:res} # dictionary for specific key and value
            resolutions.update(d) # update original resolutions dictionary

        # run the louvain clustering for the multiple resolutions
        for i in resolutions.keys():
            sc.tl.louvain(adata, resolution = resolutions[i], key_added = i)

        # Number of unique clusters generated per resolution
        res_list = list(resolutions.keys())
        [print(f'Louvain Resolution: \'{i}\' --> {adata.obs[i].nunique()} unique clusters') for i in res_list]

        if plot == True:
            sc.pl.umap(adata, color=res_list)

