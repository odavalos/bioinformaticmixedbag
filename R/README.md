# R Grabbag


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

    rig add 4.1.2-x86_64
     
M1 arm64 based (my default)

    rig default 4.1-arm64



