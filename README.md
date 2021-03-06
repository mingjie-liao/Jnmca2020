# Examples using [JuLIP](https://github.com/JuliaMolSim/JuLIP.jl)

This repository contains the numerical experiments demonstrated in JNMCA 2020. The simulation environments are `Julia@v1.3.1-v1.5.4`, `JuLIP@v0.8.13`.

## Install
`]add https://github.com/mingjie-liao/Jnmca2020.git` 

and `include` `exampleMultLJ.jl`, `exampleMultLJ3D.jl` and `plan.jl` to benchmark.

## Granular
- `multLJ.jl` is the module describing and utilizing the multi-species Lennard-Jones potential, which works as `calculator` in the numerical experiments.
- `exampleMultLJ.jl` shows the initial (randomly distributed) and stabilized configurations in two dimensions, results are saved as grad2dinit.pdf and grad2dstab.pdf.
- `exampleMultLJ3D.jl` shows the initial (randomly distributed) and stabilized configurations in three dimensions, and save the configurations in x0.dump and x1.dump respectively, which could be easily load and visualized by [ovito](https://www.ovito.org).

## Wave propogation

`plan.jl`: The wave propogation in two dimensions.
