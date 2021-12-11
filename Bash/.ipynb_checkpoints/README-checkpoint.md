## Bash

Find file names without extensions in a directory. 

`basename -s .txt *`

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


