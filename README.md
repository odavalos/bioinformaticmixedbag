# bioinformaticmixedbag


#### scRNA-seq

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





---
#### Cluster

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
#### File transfer

[ucmerced specific](https://github.com/ucmerced/merced-cluster/wiki/Transferring-Files)


FileZilla

    hostname:merced.ucmerced.edu
    username:myusername
    password:****
    port:[22](https://serverfault.com/questions/74176/what-port-does-sftp-use/167872)










