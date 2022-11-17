import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import scipy.interpolate
import numpy as np
import math 
import glob, os

import csv
import sys

import create_folder

DEM_sample_seed = [1, 2, 3] 
B = ["bridging_12_d0"]
Bridging = [0.035999999]
DX = ["dx=12d0"]
DX2 = [0.036]
C = ["XIAO", "WEAKCOUPLING"]
W = ["no_wave"]
F = ["substract_initial_force_yes", "substract_initial_force_no"]
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
nsteps = int(Nsteps / dump_every)
path_image = "Image_no_wave"
create_folder.create_folder(path_image)
Run = []
Run_wave = []

## Image folder

hh = 0.002
hh2 = 0.0025/2
err = 0.001
#Y = [DX2[0]*2, DX2[1]*4, DX2[2]*12]
color_plot = ['deepskyblue', 'steelblue', 'blue']

##  ERROR BAR PLOT - NO WAVE

DX_average_specific_section_through_time_all = []
DY_average_specific_section_through_time_all = []
DZ_average_specific_section_through_time_all = []	
DX_average_specific_section_through_time_average = []
DY_average_specific_section_through_time_average = []
DZ_average_specific_section_through_time_average = []
dx = 0
dy = 0
dz = 0

COUPLING = ["WEAKCOUPLING_substract_initial_force_yes"]
COUPLING2 = ["WEAKCOUPLING"]
COUPLING3 = ["substract_initial_force_yes"]

W = ["no_wave"]
B = ["bridging_12_d0"]
Bridging = [0.035999999]

for ii in range(len(COUPLING)):

    plot_specific = COUPLING[ii]
    coupling = COUPLING2[ii]
    substract_initial_force = COUPLING3[ii]
    for i in range(len(B)):
        bridging = B[i]
        for j in range(len(DEM_sample_seed)):
            seed_choice = DEM_sample_seed[j]
            for k in range(len(DX)):
                dx = DX[k]
                for m in range(len(W)):
                    wave = W[m]
                    if wave == "no_wave":

                        path_data = "../data/coupling/Script_coupling_no_wave/" + B[i] + "/seed_" + str(seed_choice) + "/" + dx + "/" + coupling + "/" + substract_initial_force
                        I = 1
                            
                        open_data = path_data + '/data.csv'
                        lines = []
                        with open(open_data, 'r') as file:
                            reader = csv.reader(file)
                            for line in reader:
                                lines.append(line)
                                    
                        DX_average_specific_section_through_time = []
                        DY_average_specific_section_through_time = []
                        DZ_average_specific_section_through_time = []
                            
                        for n in range(nsteps):
                            
                            DX_average_specific_section_through_time.append( float(lines[25][n]) )
                            DY_average_specific_section_through_time.append( float(lines[22][n]) )
                            DZ_average_specific_section_through_time.append( float(lines[28][n]) )
                                
                        DX_average_specific_section_through_time_all.append(DX_average_specific_section_through_time)
                        DY_average_specific_section_through_time_all.append(DY_average_specific_section_through_time)
                        DZ_average_specific_section_through_time_all.append(DZ_average_specific_section_through_time)	

    
        if plot_specific == "XIAO_substract_initial_force_yes":
            for l in range(nsteps):
                for m in range(len(DEM_sample_seed)):
                    DX_average_specific_section_through_time_all[m][l] *= 10**11
                    DY_average_specific_section_through_time_all[m][l] *= 10**11
                    DZ_average_specific_section_through_time_all[m][l] *= 10**11

        if plot_specific == "XIAO_substract_initial_force_no":
            for l in range(nsteps):
                for m in range(len(DEM_sample_seed)):
                    DX_average_specific_section_through_time_all[m][l] *= 10**6
                    DY_average_specific_section_through_time_all[m][l] *= 10**6
                    DZ_average_specific_section_through_time_all[m][l] *= 10**6
        
        if plot_specific == "WEAKCOUPLING_substract_initial_force_yes":
            for l in range(nsteps):
                for m in range(len(DEM_sample_seed)):
                    DX_average_specific_section_through_time_all[m][l] *= 10**7
                    DY_average_specific_section_through_time_all[m][l] *= 10**7
                    DZ_average_specific_section_through_time_all[m][l] *= 10**7


        if plot_specific == "WEAKCOUPLING_substract_initial_force_no":
            for l in range(nsteps):
                for m in range(len(DEM_sample_seed)):
                    DX_average_specific_section_through_time_all[m][l] *= 10**5
                    DY_average_specific_section_through_time_all[m][l] *= 10**5
                    DZ_average_specific_section_through_time_all[m][l] *= 10**5
                    
        Timestep = []
        for i in range(nsteps):
            Timestep.append(i*dump_every)
    
        # No Wave error bar
    
        fig1 = plt.figure(figsize=(3.5, 2.5))
        plt.subplot(1,1,1)
        ax1 = plt.gca()
    
        Mean_DX = []
        min_error_DX = []
        max_error_DX = []

        Mean_DY = []
        min_error_DY = []
        max_error_DY = []

        Mean_DZ = []
        min_error_DZ = []
        max_error_DZ = []

        for l in range(nsteps):
            V_DX = []
            V_DY = []
            V_DZ = []
            for m in range(len(DEM_sample_seed)):
                V_DX.append(DX_average_specific_section_through_time_all[m][l])
                V_DY.append(DY_average_specific_section_through_time_all[m][l])
                V_DZ.append(DZ_average_specific_section_through_time_all[m][l])

            min_error_DX.append(min(V_DX))
            max_error_DX.append(max(V_DX))
            Mean_DX.append(mean(V_DX))

            min_error_DY.append(min(V_DY))
            max_error_DY.append(max(V_DY))
            Mean_DY.append(mean(V_DY))
            
            min_error_DZ.append(min(V_DZ))
            max_error_DZ.append(max(V_DZ))
            Mean_DZ.append(mean(V_DZ))

        Yerr_DX = [min_error_DX, max_error_DX]
        Yerr_DY = [min_error_DY, max_error_DY]
        Yerr_DZ = [min_error_DZ, max_error_DZ]

        #ax1.errorbar(Timestep, Mean_DX, yerr=Yerr_DX, fmt='o', color='black', ecolor='gray', elinewidth=2, capsize=6, errorevery=20)
        plt.fill_between(Timestep, min_error_DX, max_error_DX, color='deepskyblue', alpha=0.2)
        ax1.plot(Timestep, Mean_DX, color='deepskyblue', linestyle = "-.", linewidth=1, label = r'$d_{x}$')
        
        #ax1.errorbar(Timestep, Mean_DY, yerr=Yerr_DY, fmt='o', color='black', ecolor='gray', elinewidth=2, capsize=6, errorevery=20)
        #print(min_error_DY, max_error_DY, Mean_DY)
        plt.fill_between(Timestep, min_error_DY, max_error_DY, color='gray', alpha=0.2)
        ax1.plot(Timestep, Mean_DY, color='slategrey', linestyle = "-", linewidth=1, label = r'$d_{y}$')
    
        #ax1.errorbar(Timestep, Mean_DZ, yerr=Yerr_DZ, fmt='o', color='black', ecolor='gray', elinewidth=2, capsize=6, errorevery=20)
        plt.fill_between(Timestep, min_error_DZ, max_error_DZ, color='blue', alpha=0.2)
        ax1.plot(Timestep, Mean_DZ, color='blue', linestyle = "--", linewidth=1, label = r'$d_{z}$')
        

        plt.xticks([0, 100000, 200000, 300000, 400000, 500000], [0, r'10$^5$', r'2.10$^5$', r'3.10$^5$', r'4.10$^5$', r'5.10$^5$' ])
    
        
        ax1.set_xlabel(r'Timestep', fontsize = 9)

        if plot_specific == "XIAO_substract_initial_force_yes":
            ax1.set_ylabel(r'Displacement ($10^{-11}$m)', fontsize = 9)
            ax1.set_title(r'Strong coupling and substraction of initial force', fontsize = 9, pad = 10)

        elif plot_specific == "XIAO_substract_initial_force_no":
            ax1.set_ylabel(r'Displacement ($10^{-6}$m)', fontsize = 9)
            ax1.set_title(r'Strong coupling and no substraction of initial force', fontsize = 9, pad = 10)

        elif plot_specific == "WEAKCOUPLING_substract_initial_force_yes":
            ax1.set_ylabel(r'Displacement ($10^{-7}$m)', fontsize = 9)
            ax1.set_title(r'Weak coupling and substraction of initial force', fontsize = 9, pad = 10)

        elif plot_specific == "WEAKCOUPLING_substract_initial_force_no":
            ax1.set_ylabel(r'Displacement ($10^{-5}$m)', fontsize = 9)
            ax1.set_title(r'Weak coupling and no substraction of initial force', fontsize = 9, pad = 10)

        else:
            ax1.set_ylabel(r'Displacement (m)', fontsize = 9)

        # plt.ylim(-2.5, 1.5)                                                                                                                                                                                  
        ax1.legend(fontsize = 8, loc='upper left')


        plt.xticks(fontsize = 8)
        plt.yticks(fontsize = 8)

        if plot_specific == "XIAO_substract_initial_force_yes":
            fig1.savefig( path_image + "/XIAO_substract_initial_force_yes.png", bbox_inches="tight")
            fig1.savefig( path_image + "/XIAO_substract_initial_force_yes.pgf", bbox_inches="tight")
            fig1.savefig( path_image + "/XIAO_substract_initial_force_yes.pdf", bbox_inches="tight")
            
        elif plot_specific == "XIAO_substract_initial_force_no":
            fig1.savefig( path_image + "/XIAO_substract_initial_force_no.png", bbox_inches="tight")
            fig1.savefig( path_image + "/XIAO_substract_initial_force_no.pgf", bbox_inches="tight")
            fig1.savefig( path_image + "/XIAO_substract_initial_force_no.pdf", bbox_inches="tight")
            
        elif plot_specific == "WEAKCOUPLING_substract_initial_force_yes":
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_yes.png", bbox_inches="tight")
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_yes.pgf", bbox_inches="tight")
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_yes.pdf", bbox_inches="tight")

        elif plot_specific == "WEAKCOUPLING_substract_initial_force_no":
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_no.png", bbox_inches="tight")
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_no.pgf", bbox_inches="tight")
            fig1.savefig( path_image + "/WEAKCOUPLING_substract_initial_force_no.pdf", bbox_inches="tight")
    
        else:
            fig1.savefig( path_image + "/No_wave.png", bbox_inches="tight")
            fig1.savefig( path_image + "/No_wave.pgf", bbox_inches="tight")
            fig1.savefig( path_image + "/No_wave.pdf", bbox_inches="tight")
    

