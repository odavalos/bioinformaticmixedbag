# bioinformaticmixedbag


## scRNA-seq

Interoperability between Seurat and Scanpy.
 - If going from Seurat to scanpy need to convert Seurat object into AnnData object.




| Seurat       | Scanpy           |
|--------------|------------------|
| nCount_RNA   | n_counts         |
| nFeature_RNA | n_genes          |
| Read10X()    | sc.read_10x_h5() |
|              |                  |
|              |                  |
|              |                  |
|              |                  |
|              |                  |
|              |                  |


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


---
## Bash

Find file names without extensions in a directory. 

`basename -s .txt *`



---
## UCM Cluster

Example sbatch job submission scripts

Example Cluster Script (Regular Job)

    # Example Cluster Script (Regular Job) (STAR Alignment)

    #!/bin/bash
    #SBATCH --mail-user=username@ucmerced.edu
    #SBATCH --mail-type=ALL
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=20
    #SBATCH --p std.q
    #SBATCH --time=0-08:45:00
    #SBATCH --output=star_blk_aln.qlog
    #SBATCH --job-name=star_blk_aln
    #SBATCH --export=ALL
    #SBATCH --oversubscribe

    cd /home/username/qsb/data/

    for R1 in *R1*
    do
       R2=${R1//R1.fastq/R2.fastq}
       read1=$R1
       read2=$R2
       base=${R1%_R1.fastq}
       _mydir=$PWD

    STAR --genomeDir /home/username/qsb/data/GCF_000001635.26_GRCm38.p6/star_index --readFilesIn $_mydir/$read1 $_mydir/$read2 --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${base}



Example Cluster Script (Array Job)

    # Example Cluster Script (Array Job) (Pulling SRA files w/parallel fastq-dump)

    #!/bin/bash
    #SBATCH --mail-user=username@ucmerced.edu
    #SBATCH --mail-type=ALL
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=20
    #SBATCH --p fast.q
    #SBATCH --time=0-3:30:00
    #SBATCH --output=GSE74148_sraarray%A_%a.qlog
    #SBATCH --job-name=GSE74148_sra
    #SBATCH --export=ALL
    #SBATCH --array=0-3
    #SBATCH --oversubscribe

    source ~/.bashrc
    cd /home/username/qsb/data/sra_files/GSE74148_Cxcr5_LCMV_HE_R # example of my dir

    sra_a=(SRR2724745 SRR2724746 SRR2724747 SRR2724748)

    echo "My SLURM_ARRAY_TASK_ID: " ${sra_a[$SLURM_ARRAY_TASK_ID]}

    prefetch -v ${sra_a[$SLURM_ARRAY_TASK_ID]}

    parallel-fastq-dump --split-files --origfmt --sra-id ${sra_a[$SLURM_ARRAY_TASK_ID]} --threads 16 --gzip








---
## File transfer

#### [ucmerced specific](https://github.com/ucmerced/merced-cluster/wiki/Transferring-Files)


FileZilla

    hostname:merced.ucmerced.edu
    username:myusername
    password:****
    port:##

[Port?](https://serverfault.com/questions/74176/what-port-does-sftp-use/167872)


---
## Conda

#### [Creating environments:](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)

Here is an example of how to create an environment with some of the essentials.


    conda create -n env_name python=3.9 pandas scipy numpy matplotlib seaborn


#### [Adding environment to jupyter kernals:](https://stackoverflow.com/a/53546634)

Example adapted from https://stackoverflow.com/a/53546634

    conda activate env_name # following the creating environment code above
    conda install ipykernel
    ipython kernel install --user --name=give_kernal_a_name
    conda deactivate
    # load up jupyter lab or jupyter notebook     


List available kernels:

`juptyer kernelspec list`



