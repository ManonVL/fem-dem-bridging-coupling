## FEM-DEM coupling

## Organization

The coupling folder contains three main folders:

- The "Script_coupling_no_wave" folder. It contains the scripts to generate the simulations presented in section 4, were no waves are generated.  In this case, three different DEM sample are considered.

- The "Script_coupling_wave" folder. It contains the scripts to generate the simulations presented in section 5, were compressive/shear waves of wavelength 0.6m and 1.2m are generated. In this case, only one DEM sample is considered.

- The folder "Plot". It contains a script which extract specific data from the simulations.

- And the "DEM" folder. It contains the script used to analyze the dem. For all the simulations, the DEM samples are given in the corresponding folder. 

## Libmultiscale simulations

To launch a libmultiscale simulation, the user must first enter in the folder of the simulation:

- cd Script_coupling_wave/bidging_12_d0/seed_1/dx=12d0/compressive_wave/wavelength_0.6/intensity_0.0005

Then, the user can launch a simulation using AMEL-sphere from libmultiscale. nsteps corresponds to the number of steps of the simulation (for the compressive wave of wavelength 0.6m nsteps = 350000 and for the shear wave of wavelength 0.6m nsteps = 500000)

- libmultiscale/build/clients/AMEL-sphere lm.config nsteps 

Results of the simulation (*.pvtu) can be vizualized with Paraview.

## Extract the results

The user can plot the results using the script main.ipynb located in the notebooks folder.

The plots are realized using the files named data.csv. This file must be generated for each simulation. 
To generate this file, the user must first unzip DEM and Plot. Then, once the simulation is over the user can create the data.csv file using:

- cd Script_coupling_wave/bidging_12_d0/seed_1/dx=12d0/compressive_wave/wavelength_0.6/intensity_0.0005
- pvpython extract_data.py

The data.csv files are included in the folders. 
