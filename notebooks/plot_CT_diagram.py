import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import scipy.interpolate
import numpy as np
import math 
import glob, os

import csv
import sys


from CT_diagram import CT_diagram
import create_folder


path_image = "Image_ct_diagram"
create_folder.create_folder(path_image)

d0 = 0.003
gouge_size = 20*d0
fe1_top = -gouge_size/2
fe2_bot = gouge_size/2

## dx = 12d0

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=12d0/compressive_wave/wavelength_0.6/intensity_5e-05/XIAO/substract_initial_force_yes"
nsteps = 350000
dump_every = 5000
wave = "compressive_wave"
name_plot = "compressive_dx=12d0_wavelength_0.6_intensity_5e-05_ct_diagram"
epsilon = 8.3e-5

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e05, epsilon, fe1_top, fe2_bot, 12*d0)

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=12d0/compressive_wave/wavelength_0.6/intensity_0.0005/XIAO/substract_initial_force_yes"
nsteps = 350000
dump_every = 5000
wave = "compressive_wave"
name_plot = "compressive_dx=12d0_wavelength_0.6_intensity_0.0005_ct_diagram"
epsilon = 8.3e-4

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e03, epsilon, fe1_top, fe2_bot, 12*d0)


path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=12d0/shear_wave/wavelength_0.6/intensity_5e-05/XIAO/substract_initial_force_yes"
nsteps = 500000
dump_every = 5000
wave = "shear_wave"
name_plot = "shear_dx=12d0_wavelength_0.6_intensity_5e-05_ct_diagram"
epsilon = 8.3e-5

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e05, epsilon, fe1_top, fe2_bot, 12*d0)

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=12d0/shear_wave/wavelength_0.6/intensity_0.0005/XIAO/substract_initial_force_yes"
nsteps = 500000
dump_every = 5000
wave = "shear_wave"
name_plot = "shear_dx=12d0_wavelength_0.6_intensity_0.0005_ct_diagram"
epsilon = 8.3e-4

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e04, epsilon, fe1_top, fe2_bot, 12*d0)

## dx = 3d0

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=3d0/compressive_wave/wavelength_0.6/intensity_5e-05/XIAO/substract_initial_force_yes"
nsteps = 350000
dump_every = 5000
wave = "compressive_wave"
name_plot = "compressive_dx=3d0_wavelength_0.6_intensity_5e-05_ct_diagram"
epsilon = 8.3e-5

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e05, epsilon, fe1_top, fe2_bot, 3*d0)

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=3d0/compressive_wave/wavelength_0.6/intensity_0.0005/XIAO/substract_initial_force_yes"
nsteps = 350000
dump_every = 5000
wave = "compressive_wave"
name_plot = "compressive_dx=3d0_wavelength_0.6_intensity_0.0005_ct_diagram"
epsilon = 8.3e-4

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e04, epsilon, fe1_top, fe2_bot, 3*d0)

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=3d0/shear_wave/wavelength_0.6/intensity_5e-05/XIAO/substract_initial_force_yes"
nsteps = 500000
dump_every = 5000
wave = "shear_wave"
name_plot = "shear_dx=3d0_wavelength_0.6_intensity_5e-05_ct_diagram"
epsilon = 8.3e-5

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e05, epsilon, fe1_top, fe2_bot, 3*d0)

path_data = "../data/coupling/Script_coupling_wave/bridging_12_d0/seed_1/dx=3d0/shear_wave/wavelength_0.6/intensity_0.0005/XIAO/substract_initial_force_yes"
nsteps = 500000
dump_every = 5000
wave = "shear_wave"
name_plot = "shear_dx=3d0_wavelength_0.6_intensity_0.0005_ct_diagram"
epsilon = 8.3e-4

CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, 1e04, epsilon, fe1_top, fe2_bot, 3*d0)

