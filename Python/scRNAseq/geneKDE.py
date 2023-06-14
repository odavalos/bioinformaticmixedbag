from sklearn.neighbors import KernelDensity
import numpy as np

def gene_kde(scan_obj, querygene, colormap = 'magma'):
    '''
    Calculate kernel density estimate and plot for a particular gene of interest.
    Function depends on KernelDensity, and numpy
    
    scan_obj: scRNAseq dataset
    querygene: Gene of interest
    colormap = cmap (colors for plotting)
    
    '''
    ##### setup #####
    # subseting scanpy object
    gene_subset = scan_obj[:, querygene]
    
    # create dataframe for subset data embeddings (UMAP)
    gene_sub_umap = pd.DataFrame(gene_subset.obsm['X_umap'], columns = ['UMAP_1', 'UMAP_2'])
    
    ##### scikit learn #####
    X = gene_subset.X.toarray() 
    kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X)
    kde_scored = kde.score_samples(X)
    
    # append log probablities to dataframe
    gene_sub_umap['log_prob'] = kde_scored * -1
    
    gene_sub_umap = gene_sub_umap.assign(COI = np.where(gene_sub_umap.log_prob >= 0, "Not COI", "COI"))
    gene_sub_umap_filt = gene_sub_umap[gene_sub_umap['COI'] == 'COI']
    
    ##### plotting #####
    
    ### single color plot ###
    
    gene_plot = plt.scatter(gene_sub_umap["UMAP_1"], gene_sub_umap["UMAP_2"],
                     c=gene_sub_umap["log_prob"], s=2, cmap=colormap) #set style options
  
    
    gene_plot = plt.title(querygene)
    gene_plot = plt.xlabel('UMAP 1')
    gene_plot = plt.ylabel('UMAP 2')
    gene_plot = plt.colorbar(label='Density')
    # gene_plot = plt.clim(0, 7)
    
    return(gene_plot)
    
