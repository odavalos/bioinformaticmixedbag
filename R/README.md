# R Grabbag


#### Colorblind Friendly Palettes

Four color blind friendly small palettes
https://thenode.biologists.com/data-visualization-with-flying-colors/research/

[Okabe_Ito](https://doi.org/10.1038/nmeth.1618)

    cbpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

[Tol Palettes](https://personal.sron.nl/~pault/)

    tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')

Colorblind palette 10 colors

    tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

Colorblind palette 9 colors

    tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')





RStudio multicursor: click and drag mouse down or up while holding down `option` 

[Modify a package function locally in R](https://stackoverflow.com/a/49277036) 

    trace("thefunction",edit=TRUE)


Use filter inside a function (Tidy Eval)

https://stackoverflow.com/a/61617180

https://stackoverflow.com/a/54948070

https://rlang.r-lib.org/reference/enquo.html


## Running multiple versions of R

**Mostly relevant if you have an M1 mac and want to run arm64 native R or an Intel based R**

- This installing and running any packages that have C++ backends
- Some bioconductor packages do not play nice with arm64

Use [Rig](https://github.com/r-lib/rig)

Here are some examples of my current installs

Intel based

    rig add 4.2.2 -a x86_64
     
M1 arm64 based (my default)

    rig default 4.1 -a arm64



## Troubleshooting 


> Using github PAT from envvar GITHUB_PAT Error: Failed to install 'unknown package' from GitHub:   
HTTP error 401.   
Bad credentials

Check this solution https://community.rstudio.com/t/error-failed-to-install-unknown-package-from-github-bad-credentials-rate-limit/113550
