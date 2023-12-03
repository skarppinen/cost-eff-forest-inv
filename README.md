# Source code and materials accompanying the paper "Planning cost-effective operational forest inventories" (2023)

This repository contains all source code and files used to draw the conclusions of the paper 
"Planning cost-effective operational forest inventories" (2023) by Santeri Karppinen, Lovisa Engberg-Sundström,
Liviu Ene and Juha Karvanen.

## Quick start 

1. Clone this repository:
```
git clone https://github.com/skarppinen/cost-eff-forest-inv.git
```

2. Install the [R programming language](https://www.r-project.org/) (code written using version 4.3.2)

- Install required R packages by executing `Rscript install-r-packages.jl` in the root folder.

3. Install the [Julia programming language](https://julialang.org/) (code written using version 1.8.5) 

- Install required Julia packages by executing `julia install-julia-packages.jl` in the root folder.

## Description of files

- `/data/test-data-4.jld2`: the dataset in the paper containing the 100 Swedish tracts.

- `/data/planning-problem-data.jld2`: the contents of `test-data-4.jld2` + the list of inventory decisions whose PoV value needs to be computed to solve the planning problem as discussed in the paper. 

- `/no-vc/output` contains the numerical results of all experiments in the paper as raw (JLD2, Julia's data format) datafiles. 

- `/src/c/random-sweep` contains the implementation of the random sweep algorithm in C. The experiments (see below) call this implementation from Julia. See separate compilation instructions below. 

- `/src/bash` a folder containing SLURM (job scheduler) files (in ZIP files) for running the experiments of the paper (on a computing cluster). These document the arguments with which the experiments of the `experiments` folder were run in order to produce the results of the paper. Note that all results are available in `no-vc/output`. 

- `/src/julia/experiments`: Folder containing the source code for launching the experiments of the paper. Each "experiment script" begins with "run-" and has a separate help menu (execute `julia run-*experiment* --help` on the command line) which lists the available arguments. 

- `/src/julia/lib`: Julia code for calling the C code, reading data, handling jobs, etc. A limited implementation of the random sweep algorithm can be found here as well. Note that the Julia code calling C requires that the associated C code has been compiled to a shared library that Julia can find. (see separate instructions for compiling below)

- `/src/julia/scripts/gen-planning-problem-data.jl` generates `/data/planning-problem-data.jld2` from `/data/test-data-4.jld2`.

- `/src/julia/scripts/gen-task-lists.jl` creates the contents of the `/src/bash/` folder.

- `/src/julia/scripts/figs` contains source code to draw the figures of the paper. Each script producing a figure begins with "plot-" and produces a figure to a newly created folder `/plots`. The remaining files are configuration for drawing the figures.

- `/figs` contains an additional figure with the tract locations shown on the map, together with the "optimal" decision obtained with the normal model.

## Compiling the C code

Running the Makefile `/src/c/random-sweep/Makefile` with `make` compiles a shared library `libpov.so` that contains the implementation of the random sweep algorithm from the paper and a function to evaluate the posterior value of inventory decision (see `pov.c` and `/src/julia/lib/pov.jl` for example of calling from Julia). By default, the GNU C compiler `gcc` is used for compilation. 
After compiling the library, move the shared library to a folder that the operating system searches when resolving shared libraries (or modify system settings such that it is found).

The C code has two dependencies:

1. [The GNU Scientific library](https://www.gnu.org/software/gsl/) (for random number generation)
2. [OpenBLAS](https://www.openblas.net/) (or some other BLAS implementation) (for linear algebra)

These dependencies have to be installed and discoverable by the operating system for compilation to succeed. 
The Makefile uses the flags `-lgsl` and `-lopenblas` to link to these libraries.

## Contact

Santeri Karppinen (skarppinen 'at' iki 'dot' fi)
