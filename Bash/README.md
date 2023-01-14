# Bash

Edit config file ~ bash shell

    emacs ~/.bashrc
    
Edit config file ~ zsh shell    

    emacs ~/.zshrc

#### Alias List

Save me the headache! Add command confirmation.

    alias rm='rm -i' # command confirmation
    alias emacs='emacs -nw'

List all the directories in current directory

    ls -ldh */

Same as above just in short version

    ls -d */
    
Return newest file in a directory ([source](https://stackoverflow.com/a/1015684))

    ls -Art | tail -n 1

Print the newest file in current directory (from https://unix.stackexchange.com/a/223305)

	cat "$(ls -rt | tail -n1)"


Find file names without extensions in a directory. 

    basename -s .txt *

List full path for all subdirectories in a directory

	find "$PWD" -type d

Get names for all files in a directory without extensions

    for i in *.ext; do basename $i .ext;done 



Renaming files using [`rename`](https://anaconda.org/bioconda/rename) (a favorite!)

	rename -v 's//GSE144236_\1/' *.csv
	# 'raw_ig_attr_Tcell.csv' renamed to 'GSE144236_raw_ig_attr_Tcell.csv'

Favorite `rename` flags

- `-n/--just-print/--dry-run`
- `-v/--verbose`


Comparing two files using `diff`
	
	diff seq1.fa seq2.fa
	-q # quickly see if files are different
	-s # shows a message
	-y # side by side differences on screen

Downloading files:

    curl -O URL
    
- Without -O curl will start dumping the downloaded file on the stdout. 
- Using -O downloads the files w/same name

Zip all files in folder in parallel

    parallel gzip ::: *


Compress a whole directory:

    tar -zcvf archive.tar.gz dir_name

	-z # compress using gzip
	-c # make archive
	-v # I like verbose settings
	-f # compressed archive file name

Extract files 

    tar -zxvf archive.tar.gz
	-x # extract files

Disk usage:

Lists sizes in readable size format for files in a directory

    du -lah

Lists sizes in readable size format in a directory

    du -h # directory/

Lists grand total for the directory

    du -sh # directory/

Get sizes of directories based on gigs

    du -ch dir/ | grep '[0-9,]G'


Get sizes of all subdirectories [(More useful examples here)](https://spapas.github.io/2018/11/12/du-disk-usage/)

    du -hs *


File Removal:

    rm file.txt # removes file

**Powerful/dangerous**; this removes all files
    
    rm *
    rm -R # removes files recursively 

Removes all files in directory ending with '.ext'

    rm *.ext

---
## Data Wrangling

Looking at a dataset in a clean manner

    head some_dataset.csv | less -S
    cat samples.txt | sed -e 's/\t/|\t|/g' | column -s, -t | less -#5 -N -S

Look at individual genes in the dataset

    less -S some_dataset.csv | grep "Tox"


Merging files with same header:

    awk 'NR==1' some_dataset.csv # this allows us to pull the header
    awk 'NR>1' data.csv | head # this prints everything but the header
    awk '(NR == 1) || (FNR > 1)' dataset1.csv dataset2.csv | (head; tail)

`(NR == 1)` includes the first line of the first file **header**, while `(FNR > 1)` skips the first line in the following files

    awk '(NR == 1) || (FNR > 1)' *.csv > merged.csv

Print stderr and stdout use `2>&1`

	kallisto quant -i /Users/username/data/mus_mus/GRCm38.p6_rna.idx -o ${base}kallisto_quant_test --threads=4 --plaintext $R1 $R2 >> kallisto.log 2>&1

Print out a log file using `2>`

	featureCounts -p -a {input.g} -t gene -g gene_id -T {threads} -o {output} {input.bams} 2> {log}

---
## Terminal Shortcuts (Mac)

Cancel the current command/line: `Ctrl+C`

Recall the deleted command: `Ctrl+Y (then Alt+Y)`

Go to beginning of the line: `Ctrl+A`

Go to end of the line: `Ctrl+E`

Remove the forward words for example, if you are middle of the command: `Ctrl+K`

Remove characters on the left, until the beginning of the word: `Ctrl+W`

To clear your entire command prompt: `Ctrl+L`

Toggle between the start of line and current cursor position: `Ctrl + XX`



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

#### My Conda setup

I like to use [Emacs text editor](https://anaconda.org/conda-forge/emacs) but need to make sure to modify ~/.bashrc

	emacs -nw ~/.bashrc
Add `alias emacs='emacs -nw'`

#### [Adding environment to jupyter kernals:](https://stackoverflow.com/a/53546634)

Example adapted from https://stackoverflow.com/a/53546634

    conda activate env_name # following the creating environment code above
    conda install ipykernel
    ipython kernel install --user --name=give_kernel_a_name
    # or use python -m ipykernel install --user --name=give_kernel_a_name
    conda deactivate
    # load up jupyter lab or jupyter notebook     


List available kernels:

`jupyter kernelspec list`

---
## Misc

IP Address Lookup

    ipconfig getifaddr en0

