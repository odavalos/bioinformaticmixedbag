## scRNA-seq

Seurat and Scanpy equivalent functions.

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




### Misc Scanpy Functions

#### Multiple louvain cluster resolution function.


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


#### Calculate kernel density estimate and plot for a particular gene of interest.

[Plot gene expression with Kernel Density Estimation (KDE) in python.](scRNAseq/gene_kde.py)

    gene_kde(adata, 'CD19', colormap = 'magma')


Example plots generated using the pbmc3k dataset from scanpy.

<img src="https://github.com/odavalos/bioinformaticmixedbag/blob/0427cb3ea2eb2932251b5bb65985cc145381941e/scRNAseq/figures/pbmc3k_kde_umap.png" width="300">

<img src="https://github.com/odavalos/bioinformaticmixedbag/blob/0427cb3ea2eb2932251b5bb65985cc145381941e/scRNAseq/figures/cd19_expression_kde_umap.png" width="200" alt="expression">

<img src="https://github.com/odavalos/bioinformaticmixedbag/blob/0427cb3ea2eb2932251b5bb65985cc145381941e/scRNAseq/figures/cd19_kde_umap.png" width="200" alt="kde">







