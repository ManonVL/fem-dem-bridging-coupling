import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import scipy.interpolate
import numpy as np
import math 
import os
import csv
import sys

from matplotlib import rc
import create_folder

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',**{'family':'serif'})
rc('text', usetex=True)

plt.rcParams.update({'figure.max_open_warning': 0})

from plot_data import  plot_data
#from pgf_delete_family import pgf_delete_family

def CT_diagram(path_data, nsteps, dump_every, wave, path_image, name_plot, I, epsilon, fe1_top, fe2_bot, mesh_size):

    # Extract data
        
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

    DX_fe1_all_steps = []
    DX_fe2_all_steps = []
    DX_dem_all_steps = []
    
    for i in range(int(nsteps/dump_every)):
        Y_fe1_all_steps.append(eval(lines[0][i]))
        Y_fe2_all_steps.append(eval(lines[1][i]))
        Y_dem_all_steps.append(eval(lines[2][i]))
        
        DY_fe1_all_steps.append(eval(lines[3][i]))
        DY_fe2_all_steps.append(eval(lines[4][i]))
        DY_dem_all_steps.append(eval(lines[5][i]))

        DX_fe1_all_steps.append(eval(lines[13][i]))
        DX_fe2_all_steps.append(eval(lines[14][i]))
        DX_dem_all_steps.append(eval(lines[15][i]))

        
    timestep = [j for j in range(0, nsteps, dump_every)]
    T1 = []
    for i in range(len(timestep)):
        T1.append([])
        for j in range(len(Y_fe1_all_steps[i])): 
            T1[i].append(timestep[i])
            
    T2 = []
    for i in range(len(timestep)):
        T2.append([])
        for j in range(len(Y_fe2_all_steps[i])): 
            T2[i].append(timestep[i])

    Tdem = []
    for i in range(len(timestep)):
        Tdem.append([])
        for j in range(len(Y_dem_all_steps[i])): 
            Tdem[i].append(timestep[i])


    for j in range(int(nsteps/dump_every)):
        for ll in range(len(DX_fe1_all_steps[j])):
            DX_fe1_all_steps[j][ll] *= I
        for ll in range(len(DY_fe1_all_steps[j])):
            DY_fe1_all_steps[j][ll] *= I
        for ll in range(len(DX_fe2_all_steps[j])):
            DX_fe2_all_steps[j][ll] *= I
        for ll in range(len(DY_fe2_all_steps[j])):
            DY_fe2_all_steps[j][ll] *= I
        for ll in range(len(DX_dem_all_steps[j])):
            DX_dem_all_steps[j][ll] *= I
        for ll in range(len(DY_dem_all_steps[j])):
            DY_dem_all_steps[j][ll] *= I

            
    Yref = []
    d0 = 0.003
    yref = fe1_top
    h = mesh_size
    while yref < fe2_bot:
        yref =  yref + h
        Yref.append(yref)

    DY_dem_all_steps_reduced = []
    DX_dem_all_steps_reduced = []
    
    for i in range(int(nsteps/dump_every)):
        
        DY_dem_all_steps_reduced.append([])
        DX_dem_all_steps_reduced.append([])
        
        for j in range(len(Yref)):
            yref = Yref[j]
            dy_reduced = 0
            dx_reduced = 0
            n = 0
            
            for k in range(len(Y_dem_all_steps[i])):
                y = Y_dem_all_steps[i][k]
                if yref - d0/2 < y < yref + d0/2:
                    #print(str(y) +" "+ str(yref))
                    dy_reduced += DY_dem_all_steps[i][k]
                    dx_reduced += DX_dem_all_steps[i][k]
                    n += 1
            if n!=0:
                dy_reduced = dy_reduced/n
                dx_reduced = dx_reduced/n
                
            DY_dem_all_steps_reduced[i].append(dy_reduced)
            DX_dem_all_steps_reduced[i].append(dx_reduced)
         
    DY_fe1 = []
    DY_fe2 = []
    DY_dem = []
    DY = []
    DX = []
    DY_fe1_all_steps_reverse = np.fliplr(np.flipud(DY_fe1_all_steps))
    DY_fe2_all_steps_reverse = np.fliplr(np.flipud(DY_fe2_all_steps))
    DY_dem_all_steps_reverse = np.fliplr(np.flipud(DY_dem_all_steps_reduced))
    DX_fe1_all_steps_reverse = np.fliplr(np.flipud(DX_fe1_all_steps))
    DX_fe2_all_steps_reverse = np.fliplr(np.flipud(DX_fe2_all_steps))
    DX_dem_all_steps_reverse = np.fliplr(np.flipud(DX_dem_all_steps_reduced))
    
    if wave == "compressive_wave":
    
        for i in range(len(timestep)):
            DY_fe1.append([])
            DY_fe2.append([])
            DY_dem.append([])
            DY.append([])
            
            for j in range(len(DY_fe2_all_steps_reverse[i])):
                DY[i].append(abs(DY_fe2_all_steps_reverse[i][j]))
                
            for j in range(len(DY_dem_all_steps_reverse[i])):
                DY[i].append(abs(DY_dem_all_steps_reverse[i][j]))
                
            for j in range(len(DY_fe1_all_steps_reverse[i])):
                DY[i].append(abs(DY_fe1_all_steps_reverse[i][j]))
            
                       
    if wave == "shear_wave":
    
        for i in range(len(timestep)):
            DX.append([])
            for j in range(len(DX_fe2_all_steps_reverse[i])):
                DX[i].append(abs(DX_fe2_all_steps_reverse[i][j]))
               
            for j in range(len(DX_dem_all_steps_reverse[i])):
                DX[i].append(abs(DX_dem_all_steps_reverse[i][j]))
                 
            for j in range(len(DX_fe1_all_steps_reverse[i])):
                DX[i].append(abs(DX_fe1_all_steps_reverse[i][j]))
            
            
    fig1 = plt.figure(figsize=(2.5, 2))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    dx, dy = 0.05, 0.05
    x = np.arange(-2.0, 2.0, dx)
    y = np.arange(-1.5, 2.0, dy)

    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    
    
    if I == 1e05:
        d = "$10^{-5}$"
    
    if I == 1e04:
        d = "$10^{-4}$"
    
    if I == 1e03:
        d = "$10^{-3}$"

    plt.title(r'$\varepsilon$ = ' + str(round(epsilon*I,1)) + '$\cdot$' + d, fontsize = 10, pad = 15)

    if wave == "compressive_wave":
        plt.xticks([-2, 0, 1.95], [-1.2, 0, 1.2])
        plt.yticks([-2, -1, 0, 1, 1.95], [0, r'10$^5$', r'2.10$^5$', r'3.10$^5$', r'4.10$^5$'])
        plt.imshow(DY, cmap='Blues_r', extent=extent, interpolation='bilinear')
    
    if wave == "shear_wave":
        plt.xticks([-2, 0, 1.95], [-1.2, 0, 1.2])
        plt.yticks([-2, -1.2, -0.4, 0.4, 1.2, 1.95], [0, r'10$^5$', r'2.10$^5$', r'3.10$^5$', r'4.10$^5$', r'5.10$^5$'])
        plt.imshow(DX, cmap='Blues_r', extent=extent, interpolation='bilinear')
    
    
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)

    ax1.set_xlabel(r'Position in the bar (m)', fontsize = 10, labelpad = 8)
    ax1.set_ylabel(r'Timestep', fontsize = 10, labelpad = 8)
    cbar = plt.colorbar()

    cbar.ax.tick_params(labelsize=8)

    if I == 1e05:
        d = "$10^{-5}$"
    
    if I == 1e04:
        d = "$10^{-4}$"
    
    if I == 1e03:
        d = "$10^{-3}$"
        
    plt.title(r'$\varepsilon$ = ' + str(round(epsilon*I,1)) + '$\cdot$' + d, fontsize = 10, pad = 15)
    
    if wave == "compressive_wave":
        cbar.set_label(r'Displacement $|u_{y}|$ (' + d + 'm)', fontsize = 10, rotation = 270, labelpad = 15)
        #cbar.set_ticks([0,1,2,3,4,5])
    if wave == "shear_wave":
        cbar.set_label(r'Displacement $|u_{x}|$ (' + d + 'm)', fontsize = 10, rotation = 270, labelpad = 15)

    path_image2 = path_image + "/PDF"
    create_folder.create_folder(path_image2)

    fig1.savefig( path_image + "/" + name_plot + ".png", bbox_inches="tight")
    fig1.savefig( path_image + "/" + name_plot + ".pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/PDF/" + name_plot + ".pdf", bbox_inches="tight")


    plt.show()
    plt.close()
