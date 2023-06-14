import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

class Abdelaal_Prepro():
    def __init__(self, min_gene, min_cell, raw_data):
        """
        A class for preprocessing scRNAseq datasets stated in \
        Abdelaal, T., Michielsen, L., Cats, D. et al. A comparison of automatic cell identification methods for single-cell RNA sequencing data. Genome Biol 20, 194 (2019). https://doi.org/10.1186/s13059-019-1795-z

        """
        self.mingene = min_gene
        self.mincell = min_cell

        if not raw_data:
            raise ValueError("No raw data provided")
        else:
            self.rawdata = raw_data
        

    def filter_genes(dataset, min_detection=0, plot:bool = False):
        """
        Filter genes with at least min_detection cells
        """
        
#         try:
#             # if it is dense
#             self.counts = np.array(self.cells.X.todense()).astype(np.float64)
#         except:
#             # if it is not dense
#             self.counts = self.cells.X
        
        sums = np.sum(dataset.X.todense() > min_detection, axis=0)
        sums = np.array(sums).flatten().reshape(-1)
        
        keep_genes = sums >= 1
        
        print(f"Remaining genes after filtering: {np.sum(keep_genes)}")
        
        filt_genes = dataset[:,keep_genes]

        if plot:

            plt.figure(figsize=(8,6))
            plt.style.use('seaborn')
            n, bins, patches = plt.hist(sums, bins = 100, facecolor='#2ab0ff', edgecolor='#e0e0e0', linewidth=0.5, alpha=0.9)

            n = n.astype('int')

            for i in range(len(patches)):
                patches[i].set_facecolor(plt.cm.RdYlBu_r(n[i]/max(n)))

            plt.xlabel('# cells in which gene is expressed')
            plt.ylabel('# of genes')
            plt.yscale('symlog')
            plt.title('Gene detection across cells')
    
        return filt_genes

    def mad_cells(dataset):
        """
        Filter out cells when the total number of detected genes 
        is bellow three MAD from the median number of detected genes per cell
        """
        try:
            # if it is not dense
            data = np.array(dataset.X.todense())
        except ValueError:
            # if it is  dense
            print("Input needs to be anndata object")
        
        total_detected_genes = np.sum(np.array(dataset.X.todense()), axis=1)


        # calculate the median number of detected genes per cell
        med = np.median(total_detected_genes)

        # calculate the median absolute deviation (mad) across all cells in the log scale
        mad = np.median(np.abs(total_detected_genes - med))
        log_mad = np.log10(mad)

        # get the absolute deviation from the median of each point 
        abs_dev = np.abs(total_detected_genes - med) / log_mad

        keep_cells = abs_dev > log_mad * 3

        print(f"Remaining cells after filtering: {np.sum(keep_cells)}")
        
        
        cells_filt = dataset[keep_cells,:]
        
        return cells_filt
      
