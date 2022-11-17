import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import scipy.interpolate
import numpy as np
import math 
import glob, os

import csv
import sys

import dem_dimensions

sys.path.insert(1, '..')
import create_folder

DEM_sample_seed = [1]
B = ["bridging_12_d0"]
Bridging = [0.035999999]
Bridging2 = [0.0361]

DX = ["dx=12d0", "dx=3d0"]
DX2 = [0.036, 0.009]
DX3= [12, 3]

C = ["XIAO"]
W = ["compressive_wave", "shear_wave"]
d0 = 0.003
L = [200*d0]
L2 = ["200"]
I = [0.0005, 0.00005]
I2 = ["5e-04", "5e-05"]
EPS = [8.3e-04, 8.3e-05]
F = ["substract_initial_force_yes"]
F2 = ["substract_initial_force ", ""]


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

path = "Image_wave_0.6"
create_folder.create_folder(path)
hh = 0.002
hh2 = 0.0025/2
err = 0.001
#Y = [DX2[0]*2, DX2[1]*4, DX2[2]*12]
color_plot_L = ['deepskyblue', 'steelblue', 'blue', 'mediumslateblue']
color_plot_I = ['deepskyblue', 'steelblue', 'blue']
#Style =  ['o', '^', '*']
Style =  ['None', 'None', 'None']
LINE = ['--','-.', ':']

## For a specific wavelength

for i in range(len(B)):

    bridging = Bridging[i]
    bridging2 = B[i]
    path1 = path + "/" + B[i]
    create_folder.create_folder(path1)
    
    for j in range(len(DEM_sample_seed)):

        seed_choice = DEM_sample_seed[j]
        path2 =  path1 + "/seed_" + str(seed_choice)
        path_dem = "../data/coupling/Script_coupling_wave/" + B[i] + "/seed_" + str(seed_choice)
        create_folder.create_folder(path2)
        
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
            
            path3 = path2 + "/" + dx
            create_folder.create_folder(path3)

            for m in range(len(W)):
                wave = W[m]
                path4 = path3 + "/" + wave
                create_folder.create_folder(path4)
                
                if wave == "compressive_wave" or wave == "shear_wave":
                    FEMLY_wave = dx_value * 12 *3
                    FEMLY = FEMLY_wave
                    for n in range(len(L)):
                        wavelength = L[n]
                        wavelength_title = L2[n]
                        path5 = path4 + "/wavelength_" + str(wavelength)
                        create_folder.create_folder(path5)

                        fig1 = plt.figure(figsize=(3.5,2.5))
                        plt.subplot(1, 1, 1)
                        ax1 = plt.gca()
                        
                        timestep = [i for i in range(0, nsteps, dump_every)]
                        
                        for o in range(len(I)):
                            intensity = I[o]
                            intensity2 = I2[o]
                            epsilon = EPS[o]
                            col = color_plot_I[o]
                            style = Style[o]
                            linest = LINE[o]
                            
                            for l in range(len(C)):
                                coupling = C[l]
                            
                                for p in range(len(F)):
                                    substract_initial_force = F[p]
                                    substract_initial_force2 = F2[p]

                                    path_data = "../data/coupling/Script_coupling_wave/" + bridging2 + "/seed_" + str(seed_choice) + "/" + dx + "/" + wave + "/wavelength_" + str(wavelength) + "/intensity_" + str(intensity) + "/" + coupling + "/" + substract_initial_force
                                    
                                    print(path_data)

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
              
                                    if wave == "shear_wave":
                                        fe1_relevant_wave_time = 85
                                        fe2_relevant_wave_time = 50
                                        open_data = path_data + '/data.csv'
                                        lines = []
                                        with open(open_data, 'r') as file:
                                            reader = csv.reader(file)
                                            for line in reader:
                                                lines.append(line)
    
                                        Y_fe1_all_steps = []
                                        Y_fe2_all_steps = []
                                        Y_dem_all_steps = []
                                        DX_fe1_all_steps = []
                                        DX_fe2_all_steps = []
                                        DX_dem_all_steps = []
    
                                        for i in range(nsteps):
                                            Y_fe1_all_steps.append(eval(lines[0][i]))
                                            Y_fe2_all_steps.append(eval(lines[1][i]))
                                            Y_dem_all_steps.append(eval(lines[2][i]))

                                        for i in range(nsteps):
                                            for j in range(len(Y_fe1_all_steps[i])):
                                                Y_fe1_all_steps[i][j] = -Y_fe1_all_steps[i][j]
                                                
                                        for i in range(nsteps):
                                            for j in range(len(Y_fe2_all_steps[i])):
                                                Y_fe2_all_steps[i][j] = -Y_fe2_all_steps[i][j]
                                                
                                        for i in range(nsteps):
                                            for j in range(len(Y_dem_all_steps[i])):
                                                Y_dem_all_steps[i][j] = -Y_dem_all_steps[i][j]

                                        for i in range(nsteps):
                                            DX_fe1_all_steps.append(eval(lines[13][i]))
                                            DX_fe2_all_steps.append(eval(lines[14][i]))
                                            DX_dem_all_steps.append(eval(lines[15][i]))

                                        for i in range(len(DX_fe1_all_steps[fe1_relevant_wave_time])):
                                            DX_fe1_all_steps[fe1_relevant_wave_time][i] = DX_fe1_all_steps[fe1_relevant_wave_time][i]/intensity
                                        for i in range(len(DX_fe2_all_steps[fe1_relevant_wave_time])):
                                            DX_fe2_all_steps[fe1_relevant_wave_time][i] = DX_fe2_all_steps[fe1_relevant_wave_time][i]/intensity
                                        for i in range(len(DX_dem_all_steps[fe1_relevant_wave_time])):    
                                            DX_dem_all_steps[fe1_relevant_wave_time][i] = DX_dem_all_steps[fe1_relevant_wave_time][i]/intensity

                                        for i in range(len(DX_fe1_all_steps[fe2_relevant_wave_time])):
                                            DX_fe1_all_steps[fe2_relevant_wave_time][i] = DX_fe1_all_steps[fe2_relevant_wave_time][i]/intensity
                                        for i in range(len(DX_fe2_all_steps[fe2_relevant_wave_time])):
                                            DX_fe2_all_steps[fe2_relevant_wave_time][i] = DX_fe2_all_steps[fe2_relevant_wave_time][i]/intensity
                                        for i in range(len(DX_dem_all_steps[fe2_relevant_wave_time])):
                                            DX_dem_all_steps[fe2_relevant_wave_time][i] = DX_dem_all_steps[fe2_relevant_wave_time][i]/intensity

                                        ax1.plot(Y_fe1_all_steps[fe1_relevant_wave_time], DX_fe1_all_steps[fe1_relevant_wave_time], color = col, linestyle = linest, marker = style, markevery = 3, linewidth=1.5, markersize=4, label = r'$\varepsilon =$ ' + "{:.1e}".format(epsilon))
                                        ax1.plot(Y_fe2_all_steps[fe1_relevant_wave_time], DX_fe2_all_steps[fe1_relevant_wave_time], color = col, linestyle = linest, marker = style,markevery = 3, linewidth=1.5, markersize=4)
                                        ax1.plot(Y_dem_all_steps[fe1_relevant_wave_time], DX_dem_all_steps[fe1_relevant_wave_time], color = col, linestyle = linest, linewidth=1.5)

                                        ax1.plot(Y_fe1_all_steps[fe2_relevant_wave_time], DX_fe1_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)
                                        ax1.plot(Y_fe2_all_steps[fe2_relevant_wave_time], DX_fe2_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)
                                        ax1.plot(Y_dem_all_steps[fe2_relevant_wave_time], DX_dem_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)


                                    if wave == "compressive_wave":
                                        fe1_relevant_wave_time = 50
                                        fe2_relevant_wave_time = 25
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

                                    
                                        ax1.plot(Y_fe1_all_steps[fe1_relevant_wave_time], DY_fe1_all_steps[fe1_relevant_wave_time], color = col,  linestyle = linest, marker = style, markevery = 3, linewidth=1.5, markersize=4, label = r'$\varepsilon =$ ' + "{:.1e}".format(epsilon))
                                        ax1.plot(Y_fe2_all_steps[fe1_relevant_wave_time], DY_fe2_all_steps[fe1_relevant_wave_time], color = col,  linestyle = linest, marker = style, markevery = 3, linewidth=1.5, markersize=4)
                                        ax1.plot(Y_dem_all_steps[fe1_relevant_wave_time], DY_dem_all_steps[fe1_relevant_wave_time], color = col,  linestyle = linest, linewidth=1.5)

                                        ax1.plot(Y_fe1_all_steps[fe2_relevant_wave_time], DY_fe1_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)
                                        ax1.plot(Y_fe2_all_steps[fe2_relevant_wave_time], DY_fe2_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)
                                        ax1.plot(Y_dem_all_steps[fe2_relevant_wave_time], DY_dem_all_steps[fe2_relevant_wave_time], color='red', linewidth=1.5)

                        ax1.set_xlabel(r'position in the bar (m)', fontsize = 9)
                        if wave == "shear_wave":
                        	ax1.set_ylabel(r'$ |u_{x} / \ u_{max}|$', fontsize = 9)
                        if wave == "compressive_wave":
                                ax1.set_ylabel(r'$ |u_{y} / \ u_{max}|$', fontsize = 9)

                        plt.title(r'Wavelength $\lambda$ =' + wavelength_title + r'$\cdot d_{0}$ and dx = ' + str(dx_value_plot) + r'$d_{0}$', fontsize = 9, pad = 10)
                        ax1.legend(fontsize = 6, loc = 'upper right')

                        plt.ylim(-0.5,1.3)
                        plt.xlim(-1.6,1.6)

                        plt.xticks(fontsize = 8)
                        plt.yticks(fontsize = 8)

                        plt.axvline(x=dem_bot,color='black',linestyle='--')
                        plt.axvline(x=dem_top,color='black',linestyle='--')

                        fig1.savefig( path5 + "/wavelength_displacement_space_" + wave + ".png", bbox_inches="tight")
                        fig1.savefig( path5 + "/wavelength_displacement_space.pgf", bbox_inches="tight")
                        fig1.savefig( path5 + "/wavelength_displacement_space.pdf", bbox_inches="tight")

                        plt.show()
                        plt.close()

