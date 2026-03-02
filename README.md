# Fredkin_Motzkin_MPS

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Fredkin_Motzkin_MPS

## Associated Publication

This repository was used to perform the numerical simulations for the paper:

Olai B. Mykland, Zhao Zhang, Exact critical exponents of the Motzkin and Fredkin Chains, arXiv: https://arxiv.org/abs/2507.14656

The data in `data/` and the figures in `plots/` correspond exactly to the results presented in the paper.

### Download
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that the data and the plots (about 1.5 MB) are included in the
   git-history and will be downloaded.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
```
which auto-activate the project and enable local path handling from DrWatson.
