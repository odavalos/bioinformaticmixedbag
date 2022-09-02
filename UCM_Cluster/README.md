# UCM Cluster

[Documentation](https://ucmerced.github.io/hpc_docs/#/)


## Data Upload and Download - Cluster

Cluster to local

    scp username@merced.ucmerced.edu:/home/username/data/fastqc/fastqc_out.html /Users/username/Downloads

Local to cluster

    scp /Users/username/data/fastqs/myfile.fastq username@merced.ucmerced.edu:/home/username/data/fastqs/myfile.fastq
    
Cluster to local (all '*.html' files in directory) 

    scp -r username@merced.ucmerced.edu:/home/username/data/fastqc/'*.html' /Users/username/data/fastqc_out

---
## Example slurm job submission scripts

Example Cluster Script (Regular Job)

    # Example Cluster Script (Regular Job) (STAR Alignment)

    #!/bin/bash
    #SBATCH --mail-user=username@ucmerced.edu
    #SBATCH --mail-type=ALL
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=20
    #SBATCH -p std.q
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

    STAR --genomeDir ~/data/GCF_000001635.26_GRCm38.p6/star_index --readFilesIn $_mydir/$read1 $_mydir/$read2 --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${base}



Example Cluster Script (Array Job)

    # Example Cluster Script (Array Job) (Pulling SRA files w/parallel fastq-dump)

    #!/bin/bash
    #SBATCH --mail-user=username@ucmerced.edu
    #SBATCH --mail-type=ALL
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=20
    #SBATCH -p fast.q
    #SBATCH --time=0-3:30:00
    #SBATCH --output=GSE74148_sraarray%A_%a.qlog
    #SBATCH --job-name=GSE74148_sra
    #SBATCH --export=ALL
    #SBATCH --array=0-3
    #SBATCH --oversubscribe

    source ~/.bashrc
    cd /home/username/data/sra_files/

    sra_a=(SRR2724745 SRR2724746 SRR2724747 SRR2724748)

    echo "My SLURM_ARRAY_TASK_ID: " ${sra_a[$SLURM_ARRAY_TASK_ID]}

    prefetch -v ${sra_a[$SLURM_ARRAY_TASK_ID]}

    parallel-fastq-dump --split-files --origfmt --sra-id ${sra_a[$SLURM_ARRAY_TASK_ID]} --threads 16 --gzip

---
## Useful Slurm commands

To delay job start add this to the job script

    #SBATCH --begin=now+1hour # this delays for 1 hour
    #SBATCH --begin=16:00 # starts at 4pm
    #SBATCH --begin=now+60 # 60 seconds (default)

If you want to test your job and find out when your job is estimated to run use (note this does not actually submit the job):

    sbatch --test-only myscript.sh

View in use and queued nodes

    squeue
    -u; e.g. flag for username to view running/queued jobs


View more detailed information about specific partitions (cleaner way than using squeue)

    sinfo
    -p, --partition; e.g. fast.q (4hr), std.q (24hr), long.q (14days)
    -t, --states; e.g. down, comp, mix, alloc (currently in use), idle (ready for use)


To cancel a job

    scancel jobid

To cancel all the jobs for a user

    scancel -u username


