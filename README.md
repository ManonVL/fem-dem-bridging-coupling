# FEM-DEM coupling 

## Getting started

In this repository you will find the scripts and data referring to the paper "FEM-DEM bridging coupling for the modeling of gouge" (put doi once available).

This repository is composed of two main folders:

- The data folder, which is composed of:

	- The elastic_properties folder, in which you will find the scripts and the data used to obtain the results presented in the section 3 of the paper.
	
	- The coupling folder, in which you will find the scripts and the data used to obtain the results presented in the sections 4 and 5 of the paper.

		Each sub folder contains its own README.md file.

- The notebook folder, in which you will find the main.ipynb file which allow the user to reproduce the plots of the paper.

## Tools

The user will need python, lammps, paraview, and libmultiscale to properly launch all the scripts. 

The version of libmultiscale used in the paper to launch the simulations is : fde0d490b83e9ff4c613322841408b733affa8a2
