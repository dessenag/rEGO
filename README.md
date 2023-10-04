# rEGO [![DOI](https://zenodo.org/badge/699360255.svg)](https://zenodo.org/badge/latestdoi/699360255)
A [tutorial](/Tutorial) for the refined Efficient Optimisation algorithm


This repo contains the code for the article:

G. Dessena, D. I. Ignatyev, J. F. Whidborne, and L. Zanotti Fragonara, ‘A global-local metamodeling technique for model updating’, [Accepted] Computer Methods in Applied Mechanics and Engineering, 2023.

When using this code for your work or research please cite the following:

1.	G. Dessena, D. I. Ignatyev, J. F. Whidborne, and L. Zanotti Fragonara, ‘A global-local metamodeling technique for model updating’, [Accepted] Computer Methods in Applied Mechanics and Engineering, 2023.
2.	G. Dessena, D. I. Ignatyev, J. F. Whidborne, and L. Zanotti Fragonara, ‘A Kriging Approach to Model Updating for Damage Detection’, EWSHM 2022, Lecture Notes in Civil Engineering. Springer International Publishing, pp. 245–255, Jun. 16, 2022. doi: [10.1007/978-3-031-07258-1_26](https://doi.org/10.1007/978-3-031-07258-1_26).
3.	G. Dessena, rEGO – A tutorial on the refined Efficient Global Optimisation, GitHub, Oct. 4, 2023. doi: [10.5281/zenodo.8406031](https://doi.org/10.5281/zenodo.8406031)


The rEGO is a refined version of the well-known EGO and allows for a better use of computational power for single-objective problems. The repository contains also the model, built in [1] for the three-storey structure of the EI at LANL. 
The rEGO workflow is a novelty and the code is built over two existing platforms which are forked under [Utilities](/Utilities):

- DACE Toolbox, Hans Bruun Nielsen, Søren Nymand and Lophaven Jacob Søndergaard - https://github.com/psbiomech/dace-toolbox-source
- Single_objective_EGO_algorithms, Qi Zhang - https://github.com/109902249/Single_objective_EGO_algorithms

## Table of Contents
- [rEGO.m](/rEGO.m): main file of the rEGO algorithm
- [Tutorial](/Tutorial): MATLAB tutorials in .m, .mlx, and pdf format for rEGO on two test functions
- [Utilities](/Utilities): forks, test functions, and other misc functions.
- [LANL_3SS](/LANL_3SS): a numerical model of the Three Storey Structure from the Engineering Institute at LANL (https://doi.org/10.2172/961604)

## Repository Citation

G. Dessena, rEGO – A tutorial on the refined Efficient Global Optimisation, GitHub, Oct. 4, 2023. doi: [10.5281/zenodo.8406031](https://doi.org/10.5281/zenodo.8406031)