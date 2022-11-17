import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import scipy.interpolate
import numpy as np
import math 
import glob, os

import csv
import sys

sys.path.insert(1, '../DEM')
import dem_dimensions

sys.path.insert(1, '..')
import create_folder

from matplotlib import rc

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',**{'family':'serif'})
rc('text', usetex=True)

sys.path.append("Generate_plots_and_videos")

plt.rcParams.update({'figure.max_open_warning': 0})

DEM_sample_seed = [1]
B = ["bridging_12_d0"]
Bridging = [0.035999999]
Bridging2 = [0.0361]

#DX = ["dx=10d0", "dx=3d0"]
#DX2 = [0.036, 0.009]
#DX3= [12, 3]

DX = ["dx=12d0"]
DX2 = [0.036]
DX3= [12]

C = ["XIAO"]
W = ["compressive_wave"]
d0 = 0.003
#L = [400*d0, 200*d0, 100*d0, 50*d0, 30*d0]
#L2 = [400, 200, 100, 50, 30]
L = [400*d0, 200*d0]
L2 = [400, 200]
intensity = 0.00005
I2 = ["5e-05"]
EPS = [4.2e-05, 8.3e-05, 1.7e-04]

substract_initial_force =  "substract_initial_force_yes"
coupling = "XIAO" 
#color_plot_I = ['deepskyblue', 'steelblue']
#color_plot_I2 = ['red', 'orange']

color_plot_I = ['coral', 'deepskyblue']
color_plot_I2 = ['red', 'navy']

# color_plot_I = ['coral', 'steelblue']
# color_plot_I2 = ['red', 'navy']


#Style =  ['o', '^', '*']
Style =  ['None', 'None']
#LINE = ['-','--']

# Geometry parameters

d0 = 0.003
dump_every = 5000
FEMLY_wave = 1.2
FEMLY_no_wave = 0.4
gouge_size = 20*d0

# FEM material properties

name_fem = "dem_composed_of_glass_bead"
rho_fem = 1643  # density in kg/m^3
E_fem = 977000000  # young's in Pa
nu_fem = 0.255
eigen_grad_u = [[0.0025, 0, 0], [0, 0.0025, 0], [0, 0, 0.0025]]

# DEM

E = 50*10**9
nu = 0.3
friction_between_grains = 0.5
coeff_of_restitution = 0.9
density = 2500
dmin = 0.0025
dmax = 0.0035
dt = 1*10**-8
d1 = 0.0025
dmin = 0.0025
dmax = 0.0035
d0 = 0.003
d2 = 0.0035
de = 0.0035
gamma_n = 0.2
gamma_t = 'NULL'
damping = 0.05
kn = (2/3)*E/(1-nu**2)
kt = 2*E/((1+nu)*(2-nu))
xmu = friction_between_grains
dampflag = 1
pressure = 5e6
visc = 0.5

Nsteps = 500000
dump_every = 5000
nsteps = int( Nsteps / dump_every )
Run = []
Run_wave = []

## Image folder

path = "Image_several_wavelength"
create_folder.create_folder(path)
hh = 0.002
hh2 = 0.0025/2
err = 0.001
#Y = [DX2[0]*2, DX2[1]*4, DX2[2]*12]

## For a specific wavelength

for i in range(len(B)):

    bridging = Bridging[i]
    bridging2 = B[i]
   
    for j in range(len(DEM_sample_seed)):

        seed_choice = DEM_sample_seed[j]
        path_dem = "../data/coupling/Script_coupling_wave/" + B[i] + "/seed_" + str(seed_choice)
        
        DEMLX, DEMLY, DEMLZ = dem_dimensions.dem_dimensions(path_dem, "confined_5.0_MPa_duplicated_position.data")
        DEMLX = round(DEMLX,6)
        DEMLY = round(DEMLY,6)
        DEMLZ = round(DEMLZ,6)
        FEMLX = DEMLX
        FEMLZ = DEMLZ

        for k in range(len(DX)):
            dx = DX[k]
            dx_value = DX2[k]
            dx_value_plot = DX3[k]
            Nelement_x = round(FEMLX/dx_value) + 1

            for m in range(len(W)):
                wave = W[m]

                fig1 = plt.figure(figsize=(3.5,2.5))
                plt.subplot(1, 1, 1)
                ax1 = plt.gca()

                timestep = [i for i in range(0, nsteps, dump_every)]
                
                if wave == "compressive_wave":
                    FEMLY_wave = dx_value * 12 *3
                    FEMLY = FEMLY_wave
                    for n in range(len(L)):
                        col = color_plot_I[n]
                        col2 = color_plot_I2[n]
                        style = Style[n]
                        #linest = LINE[n]
                        epsilon = EPS[n] 
                        wavelength = L[n]
                        wavelength_title = L2[n]
                        
                        path_data = "../data/coupling/Script_coupling_wave/" + bridging2 + "/seed_" + str(seed_choice) + "/" + dx + "/" + wave + "/wavelength_" + str(wavelength) + "/intensity_" + str(intensity) + "/" + coupling + "/" + substract_initial_force
                       

                        fe1_bot = -FEMLY - gouge_size/2
                        fe1_top = -gouge_size/2
                        fe2_bot = gouge_size/2
                        fe2_top = FEMLY + gouge_size/2
                        dem_bot = -DEMLY/2
                        dem_top = DEMLY/2
                        y_value_fem = fe2_bot + dx_value
                        hh = 0.002
                        hh2 = 0.0025/2
                        err = 0.001

                        if wave == "compressive_wave" and wavelength == 400*d0:
                            fe1_relevant_wave_time = 100
                            fe2_relevant_wave_time = 50
                            nsteps = 130
                        if wave == "compressive_wave" and wavelength != 400*d0:
                            fe1_relevant_wave_time = 50
                            fe2_relevant_wave_time = 25
                            nsteps = 100
                        open_data = path_data + '/data.csv'
                        lines = []
                        with open(open_data, 'r') as file:
                            reader = csv.reader(file)
                            for line in reader:
                                lines.append(line)
                                                
                        Y_fe1_all_steps = []
                        Y_fe2_all_steps = []
                        Y_dem_all_steps = []
                        DY_fe1_all_steps = []
                        DY_fe2_all_steps = []
                        DY_dem_all_steps = []
    
                        for i in range(nsteps):
                            Y_fe1_all_steps.append(eval(lines[0][i]))
                            Y_fe2_all_steps.append(eval(lines[1][i]))
                            Y_dem_all_steps.append(eval(lines[2][i]))
                                            
                        for i in range(nsteps):
                            DY_fe1_all_steps.append(eval(lines[3][i]))
                            DY_fe2_all_steps.append(eval(lines[4][i]))
                            DY_dem_all_steps.append(eval(lines[5][i]))

                        for i in range(nsteps):
                            for j in range(len(Y_fe1_all_steps[i])):
                                Y_fe1_all_steps[i][j] = -Y_fe1_all_steps[i][j]
                                                
                        for i in range(nsteps):
                            for j in range(len(Y_fe2_all_steps[i])):
                                Y_fe2_all_steps[i][j] = -Y_fe2_all_steps[i][j]
                                                
                        for i in range(nsteps):
                            for j in range(len(Y_dem_all_steps[i])):
                                Y_dem_all_steps[i][j] = -Y_dem_all_steps[i][j]

                        for i in range(len(DY_fe1_all_steps[fe1_relevant_wave_time])):
                            DY_fe1_all_steps[fe1_relevant_wave_time][i] = -DY_fe1_all_steps[fe1_relevant_wave_time][i]/intensity
                        for i in range(len(DY_fe2_all_steps[fe1_relevant_wave_time])):
                            DY_fe2_all_steps[fe1_relevant_wave_time][i] = -DY_fe2_all_steps[fe1_relevant_wave_time][i]/intensity
                        for i in range(len(DY_dem_all_steps[fe1_relevant_wave_time])):    
                            DY_dem_all_steps[fe1_relevant_wave_time][i] = -DY_dem_all_steps[fe1_relevant_wave_time][i]/intensity
                                        
                        for i in range(len(DY_fe1_all_steps[fe2_relevant_wave_time])):
                            DY_fe1_all_steps[fe2_relevant_wave_time][i] = -DY_fe1_all_steps[fe2_relevant_wave_time][i]/intensity
                        for i in range(len(DY_fe2_all_steps[fe2_relevant_wave_time])):
                            DY_fe2_all_steps[fe2_relevant_wave_time][i] = -DY_fe2_all_steps[fe2_relevant_wave_time][i]/intensity
                        for i in range(len(DY_dem_all_steps[fe2_relevant_wave_time])):
                            DY_dem_all_steps[fe2_relevant_wave_time][i] = -DY_dem_all_steps[fe2_relevant_wave_time][i]/intensity

                        if wavelength == 400*d0:
                            line1, = ax1.plot(Y_fe1_all_steps[fe2_relevant_wave_time], DY_fe1_all_steps[fe2_relevant_wave_time], color = col2, linestyle = '-', linewidth=1.5, label = r'$\lambda =$ ' + str(wavelength_title) + r'$d_{0}$')
                        if wavelength == 200*d0:
                            line11, = ax1.plot(Y_fe1_all_steps[fe2_relevant_wave_time], DY_fe1_all_steps[fe2_relevant_wave_time], color = col2, linestyle = '-', linewidth=1.5, label = r'$\lambda =$ ' + str(wavelength_title) + r'$d_{0}$')
                        ax1.plot(Y_fe2_all_steps[fe2_relevant_wave_time], DY_fe2_all_steps[fe2_relevant_wave_time], color = col2, linestyle = '-', linewidth=1.5)
                        ax1.plot(Y_dem_all_steps[fe2_relevant_wave_time], DY_dem_all_steps[fe2_relevant_wave_time], color = col2, linestyle = '-', linewidth=1.5)
            
                        if wavelength == 400*d0:
                            line2, = ax1.plot(Y_fe1_all_steps[fe1_relevant_wave_time], DY_fe1_all_steps[fe1_relevant_wave_time], color = col, linestyle = '--', marker = style, markevery = 3, linewidth=1.5, markersize=4, label = r'$\lambda =$ ' + str(wavelength_title) + r'$d_{0}$')
                        if wavelength == 200*d0:  
                            line22, = ax1.plot(Y_fe1_all_steps[fe1_relevant_wave_time], DY_fe1_all_steps[fe1_relevant_wave_time], color = col, linestyle = '--', marker = style, markevery = 3, linewidth=1.5, markersize=4, label = r'$\lambda =$ ' + str(wavelength_title) + r'$d_{0}$')
                        ax1.plot(Y_fe2_all_steps[fe1_relevant_wave_time], DY_fe2_all_steps[fe1_relevant_wave_time], color = col, linestyle = '--', marker = style, markevery = 3, linewidth=1.5, markersize=4)
                        ax1.plot(Y_dem_all_steps[fe1_relevant_wave_time], DY_dem_all_steps[fe1_relevant_wave_time], color = col, linestyle = '--', linewidth=1.5)

                        
                ax1.set_xlabel(r'position in the bar (m)', fontsize = 9)
                if wave == "shear_wave":
                    ax1.set_ylabel(r'$ |u_{x} / \ u_{x,\mathrm{max}}|$', fontsize = 9)
                if wave == "compressive_wave":
                    ax1.set_ylabel(r'$ |u_{y} / \ u_{y,\mathrm{max}}|$', fontsize = 9)

                plt.title(r'Intensity $I$ =' + " {:.1e}".format(intensity) + r"(m) ", fontsize = 9, pad = 10)
      
                first_legend = ax1.legend(handles=[line1, line11], fontsize = 6, loc='upper left', title = 'Before:')
                #first_legend.get_title().set_fontsize('6')
                first_legend = ax1.legend(handles=[line1, line11], fontsize = 6, loc='upper left')
                ax1.add_artist(first_legend)
                #second_legend = ax1.legend(handles=[line2, line22], fontsize = 6, loc='upper right', title = 'After:')
                #second_legend.get_title().set_fontsize('6')
                second_legend = ax1.legend(handles=[line2, line22], fontsize = 6, loc='upper right')
                
                #plt.ylim(-0.5,1.3)
                plt.xlim(-3.2,3.2)
                plt.xticks([-3, -2, -1, 0, 1, 2, 3], [-3, -2, -1, 0, 1, 2, 3])
                plt.xticks(fontsize = 8)
                plt.yticks(fontsize = 8)

                plt.axvline(x=dem_bot,color='black',linestyle='--', linewidth=1.5)
                plt.axvline(x=dem_top,color='black',linestyle='--', linewidth=1.5)

                fig1.savefig( path + "/compressive_new.png", bbox_inches="tight")
                fig1.savefig( path + "/compressive_new.pgf", bbox_inches="tight")
                fig1.savefig( path + "/compressive_new.pdf", bbox_inches="tight")
                #pgf_delete_family( path_image + "/displacement_space.pgf", path_image + "/displacement_space_latex.pgf")
                plt.show()
                plt.close()


