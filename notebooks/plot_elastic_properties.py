import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pylab import *
import numpy as np
import math 
import os
import csv

from import_atoms import import_atoms
from density_tot import density_tot
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',**{'family':'serif'})
rc('text', usetex=True)

plt.rcParams.update({'figure.max_open_warning': 0})

P = [5e6]
d0 = 0.003
dmean = 0.003
Size = [2*d0, 3*d0, 6*d0, 9*d0, 12*d0, 15*d0, 18*d0, 21*d0, 24*d0, 27*d0]
Size_title = ["2", "3", "6", "9", "12", "15", "18", "21", "24", "27"]
Seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19]

path_image = "Image_elastic_properties"
Table = []
E_glass_bead = 50*10**9

try:
    os.mkdir(path_image)
except OSError:
    print ("directory already exists")
else:
    print ("Successfully created the directory %s" % path_image)


for i in range(len(P)):

    pressure = P[i]
    
   
    Bulk_modulus_mean = []
    Shear_modulus1_mean = []
    Shear_modulus2_mean = []
    Shear_modulus_average_mean = []
    Poisson_mean = []
    Anisotropy_coefficient_mean = []
    Young_modulus_mean = []
    Size_box_mean = []
    Density_mean = []
    
    Bulk_modulus_all = []
    Shear_modulus1_all = []
    Shear_modulus2_all = []
    Shear_modulus_average_all = []
    Poisson_all = []
    Anisotropy_coefficient_all = []
    Young_modulus_all = []
    Size_box_all = []
    Density_all = []
    
    Bulk_modulus_error_bar = []
    Shear_modulus1_error_bar = []
    Shear_modulus2_error_bar = []
    Shear_modulus_average_error_bar = []
    Poisson_error_bar = []
    Anisotropy_coefficient_error_bar = []
    Young_modulus_error_bar = []
    Density_error_bar = []
    Length_box_error_bar = []

    n = -1
    
    for j in range(len(Size)):
        
        Lx = Size[j]    
        s_title = Size_title[j]
        
        C11_average = 0
        C12_average = 0
        C13_average = 0
        C14_average = 0
        C15_average = 0
        C16_average = 0
        C26_average = 0
        C21_average = 0
        C22_average = 0
        C23_average = 0
        C24_average = 0
        C25_average = 0
        C26_average = 0
        C31_average = 0
        C32_average = 0
        C33_average = 0
        C34_average = 0
        C35_average = 0
        C36_average = 0
        C41_average = 0
        C42_average = 0
        C43_average = 0
        C44_average = 0
        C45_average = 0
        C46_average = 0
        C51_average = 0
        C52_average = 0
        C53_average = 0
        C54_average = 0
        C55_average = 0
        C56_average = 0
        C61_average = 0
        C62_average = 0
        C63_average = 0
        C64_average = 0
        C65_average = 0
        C66_average = 0
        
        bulk_modulus_average = 0
        shear_modulus1_average = 0
        shear_modulus2_average = 0
        poisson_ratio_average = 0
        anisotropy_coefficient_average = 0
        young_modulus_average = 0
        density_average = 0
        L_average = 0   
        
        Bulk_modulus_error_bar.append([])
        Shear_modulus1_error_bar.append([])
        Shear_modulus2_error_bar.append([])
        Shear_modulus_average_error_bar.append([])
        Poisson_error_bar.append([])
        Anisotropy_coefficient_error_bar.append([])
        Young_modulus_error_bar.append([])
        Density_error_bar.append([])
        Length_box_error_bar.append([])
    
        Bulk_modulus_mean.append([])
        Shear_modulus1_mean.append([])
        Shear_modulus2_mean.append([])
        Shear_modulus_average_mean.append([])
        Poisson_mean.append([])
        Anisotropy_coefficient_mean.append([])
        Young_modulus_mean.append([])
        Size_box_mean.append([])
        Density_mean.append([])
    
        Bulk_modulus_specific_size = []
        Shear_modulus1_specific_size = []
        Shear_modulus2_specific_size = []
        Shear_modulus_average_specific_size = []
        Poisson_specific_size = []
        Anisotropy_coefficient_specific_size = []
        Young_modulus_specific_size = []
        Density_specific_size = []
        L_box = []
        
        n = n + 1
        
        for k in range(len(Seed)):

            seedChoice = Seed[k]
                  
            path = "../data/elastic_properties/box_" + str(Lx) + "/seed_" + str(seedChoice) + "/pressure_" + str(pressure*10**-6) + "/elastic_properties"

            # Elastic properties
            open1 = "../data/elastic_properties/box_" + str(Lx) + "/seed_" + str(seedChoice) + "/pressure_" + str(pressure*10**-6) + "/elastic_properties/stiffness_matrix_sym.txt"
            f1 = open(open1, "r")
            A1 = []
            A1 = [ line.split() for line in f1]
            del A1[0:1]
           
            
            open2 = "../data/elastic_properties/box_" + str(Lx) + "/seed_" + str(seedChoice) + "/pressure_" + str(pressure*10**-6) + "/elastic_properties/properties.txt"
            f2 = open(open2, "r")
            A2 = []
            A2 = [ line.split() for line in f2]
            del A2[0:1]
            
            bulk_modulus_average = bulk_modulus_average + float(A2[0][0])
            shear_modulus1_average = shear_modulus1_average + float(A2[0][1])
            shear_modulus2_average = shear_modulus2_average + float(A2[0][2])
            poisson_ratio_average = poisson_ratio_average + float(A2[0][3])
            anisotropy_coefficient_average = anisotropy_coefficient_average + float(A2[0][4])
            young_modulus_average = young_modulus_average + float(A2[0][5])
            
            Bulk_modulus_specific_size.append(float(A2[0][0]))
            Shear_modulus1_specific_size.append(float(A2[0][1]))
            Shear_modulus2_specific_size.append(float(A2[0][2]))
            Shear_modulus_average_specific_size.append((float(A2[0][1])+float(A2[0][2]))/2)
            Poisson_specific_size.append(float(A2[0][3]))
            Anisotropy_coefficient_specific_size.append(float(A2[0][4]))
            Young_modulus_specific_size.append(float(A2[0][5]))

            # Box_size
            input_file = "../data/elastic_properties/box_" + str(Lx) + "/seed_" + str(seedChoice) + "/pressure_" + str(pressure*10**-6) + "/confined_" + str(pressure*10**-6) + "_MPa.data"
            atoms, boundaries = import_atoms(input_file)
            L = boundaries[1] - boundaries[0]
            L_box.append(L)
            L_average = L_average + L
            
            # Total density 
            d_tot = density_tot(input_file)
            Density_specific_size.append(d_tot)
            density_average = density_average + d_tot

            
            ## Stiffness matrix sym
            C11_average = C11_average + float(A1[0][0])
            C12_average = C12_average + float(A1[0][1])
            C13_average = C13_average + float(A1[0][2])
            C14_average = C14_average + float(A1[0][3])
            C15_average = C15_average + float(A1[0][4])
            C16_average = C16_average + float(A1[0][5])
            
            C21_average = C21_average + float(A1[0][6])
            C22_average = C22_average + float(A1[0][7])
            C23_average = C23_average + float(A1[0][8])
            C24_average = C24_average + float(A1[0][9])
            C25_average = C25_average + float(A1[0][10])
            C26_average = C26_average + float(A1[0][11])
            
            C31_average = C31_average + float(A1[0][12])
            C32_average = C32_average + float(A1[0][13])
            C33_average = C33_average + float(A1[0][14])
            C34_average = C34_average + float(A1[0][15])
            C35_average = C35_average + float(A1[0][16])
            C36_average = C36_average + float(A1[0][17])
            
            C41_average = C41_average + float(A1[0][18])
            C42_average = C42_average + float(A1[0][19])
            C43_average = C43_average + float(A1[0][20])
            C44_average = C44_average + float(A1[0][21])
            C45_average = C45_average + float(A1[0][22])
            C46_average = C46_average + float(A1[0][23])
            
            C51_average = C51_average + float(A1[0][24])
            C52_average = C52_average + float(A1[0][25])
            C53_average = C53_average + float(A1[0][26])
            C54_average = C54_average + float(A1[0][27])
            C55_average = C55_average + float(A1[0][28])
            C56_average = C56_average + float(A1[0][29])
            
            C61_average = C61_average + float(A1[0][30])
            C62_average = C62_average + float(A1[0][31])
            C63_average = C63_average + float(A1[0][32])
            C64_average = C64_average + float(A1[0][33])
            C65_average = C65_average + float(A1[0][34])
            C66_average = C66_average + float(A1[0][35])           
            
            
            C = []
            C.append( [float(A1[0][0]), float(A1[0][1]), float(A1[0][2]), float(A1[0][3]), float(A1[0][4]), float(A1[0][5])] )
            C.append( [float(A1[0][6]), float(A1[0][7]), float(A1[0][8]), float(A1[0][9]), float(A1[0][10]), float(A1[0][11])] )
            C.append( [float(A1[0][12]), float(A1[0][13]), float(A1[0][14]), float(A1[0][15]), float(A1[0][16]), float(A1[0][17])] )
            C.append( [float(A1[0][18]), float(A1[0][19]), float(A1[0][20]), float(A1[0][21]), float(A1[0][22]), float(A1[0][23])] )
            C.append( [float(A1[0][24]), float(A1[0][25]), float(A1[0][26]), float(A1[0][27]), float(A1[0][28]), float(A1[0][29])] )
            C.append( [float(A1[0][30]), float(A1[0][31]), float(A1[0][32]), float(A1[0][33]), float(A1[0][34]), float(A1[0][35])] )

            # Change units of the stiffness marix Pa --> GPa
            for i in range(6):
               for j in range(6):
                  C[i][j] = C[i][j]*10**-9
 
        
        Bulk_modulus_all.append(Bulk_modulus_specific_size)
        Shear_modulus1_all.append(Shear_modulus1_specific_size)
        Shear_modulus2_all.append(Shear_modulus2_specific_size)
        Shear_modulus_average_all.append(Shear_modulus_average_specific_size)
        Poisson_all.append(Poisson_specific_size)
        Anisotropy_coefficient_all.append(Anisotropy_coefficient_specific_size)
        Young_modulus_all.append(Young_modulus_specific_size)
        Density_all.append(Density_specific_size)
        Size_box_all.append(L_box)
        
        bulk_modulus_average = bulk_modulus_average/len(Seed)
        shear_modulus1_average = shear_modulus1_average/len(Seed)
        shear_modulus2_average = shear_modulus2_average/len(Seed)
        poisson_ratio_average = poisson_ratio_average/len(Seed)
        anisotropy_coefficient_average = anisotropy_coefficient_average/len(Seed)
        young_modulus_average = young_modulus_average/len(Seed)
        density_average = density_average/len(Seed)
        L_average = L_average/len(Seed)
        
        Bulk_modulus_mean[n].append(bulk_modulus_average)
        Shear_modulus1_mean[n].append(shear_modulus1_average)
        Shear_modulus2_mean[n].append(shear_modulus2_average)
        Shear_modulus_average_mean[n].append((shear_modulus1_average+shear_modulus2_average)/2)
        Poisson_mean[n].append(poisson_ratio_average)
        Anisotropy_coefficient_mean[n].append(anisotropy_coefficient_average)
        Young_modulus_mean[n].append(young_modulus_average)
        Size_box_mean[n].append(L_average)
        Density_mean[n].append(density_average)
        bulk_error = 0
        shear1_error = 0
        shear2_error = 0
        poisson_error = 0
        anisotropy_error = 0
        young_error = 0
        density_error = 0
        length_box_error = 0
        for l in range(len(Seed)):
            bulk_error = bulk_error + (Bulk_modulus_specific_size[l] - bulk_modulus_average)**2
            shear1_error = shear1_error + (Shear_modulus1_specific_size[l] - shear_modulus1_average )**2
            shear2_error = shear2_error + (Shear_modulus2_specific_size[l] - shear_modulus2_average)**2
            poisson_error = poisson_error + (Poisson_specific_size[l] - poisson_ratio_average)**2
            anisotropy_error = anisotropy_error + (Anisotropy_coefficient_specific_size[l] - anisotropy_coefficient_average)**2
            young_error = young_error + (Young_modulus_specific_size[l] - young_modulus_average)**2
            density_error = density_error + (Density_specific_size[l] - density_average)**2
            length_box_error = length_box_error + (L_box[l]-L_average)**2
            
        bulk_error = math.sqrt(bulk_error/len(Bulk_modulus_specific_size))
        shear1_error = math.sqrt(shear1_error/len(Shear_modulus1_specific_size))
        shear2_error = math.sqrt(shear2_error/len(Shear_modulus2_specific_size))
        poisson_error = math.sqrt(poisson_error/len(Poisson_specific_size))
        young_error = math.sqrt(young_error/len(Young_modulus_specific_size))
        density_error = math.sqrt(density_error/len(Density_specific_size))
        length_box_error = math.sqrt(length_box_error/len(L_box))
        
        Bulk_modulus_error_bar[n].append(bulk_error)
        Shear_modulus1_error_bar[n].append(shear1_error)
        Shear_modulus2_error_bar[n].append(shear2_error)
        Shear_modulus_average_error_bar[n].append((shear1_error+shear2_error)/2)
        Poisson_error_bar[n].append(poisson_error)
        Anisotropy_coefficient_error_bar[n].append(anisotropy_error)
        Young_modulus_error_bar[n].append(young_error)
        Density_error_bar[n].append(density_error)
        Length_box_error_bar[n].append(length_box_error)
        
        C11_average = C11_average/len(Seed)
        C12_average = C12_average/len(Seed)
        C13_average = C13_average/len(Seed)
        C14_average = C14_average/len(Seed)
        C15_average = C15_average/len(Seed)
        C16_average = C16_average/len(Seed)
        
        C21_average = C21_average/len(Seed)
        C22_average = C22_average/len(Seed)
        C23_average = C23_average/len(Seed)
        C24_average = C24_average/len(Seed)
        C25_average = C25_average/len(Seed)
        C26_average = C26_average/len(Seed)
        
        C31_average = C31_average/len(Seed)
        C32_average = C32_average/len(Seed)
        C33_average = C33_average/len(Seed)
        C34_average = C34_average/len(Seed)
        C35_average = C35_average/len(Seed)
        C36_average = C36_average/len(Seed)
        
        C41_average = C41_average/len(Seed)
        C42_average = C42_average/len(Seed)
        C43_average = C43_average/len(Seed)
        C44_average = C44_average/len(Seed)
        C45_average = C45_average/len(Seed)
        C46_average = C46_average/len(Seed)
        
        C51_average = C51_average/len(Seed)
        C52_average = C52_average/len(Seed)
        C53_average = C53_average/len(Seed)
        C54_average = C54_average/len(Seed)
        C55_average = C55_average/len(Seed)
        C56_average = C56_average/len(Seed)
        
        C61_average = C61_average/len(Seed)
        C62_average = C62_average/len(Seed)
        C63_average = C63_average/len(Seed)
        C64_average = C64_average/len(Seed)
        C65_average = C65_average/len(Seed)
        C66_average = C66_average/len(Seed)
        
        C = []
        C.append( [C11_average, C12_average, C13_average, C14_average, C15_average, C16_average] )
        C.append( [C21_average, C22_average, C23_average, C24_average, C25_average, C26_average] )
        C.append( [C31_average, C32_average, C33_average, C34_average, C35_average, C36_average] )
        C.append( [C41_average, C42_average, C43_average, C44_average, C45_average, C46_average] )
        C.append( [C51_average, C52_average, C53_average, C54_average, C55_average, C56_average] )
        C.append( [C61_average, C62_average, C63_average, C64_average, C65_average, C66_average] )

        # Change units of the stiffness marix Pa --> GPa
        for i in range(6):
           for j in range(6):
               C[i][j] = C[i][j]*10**-9

        Table.append([C11_average, C12_average, C13_average, C14_average, C15_average, C16_average, C21_average, C22_average, C23_average, C24_average, C25_average, C26_average, C31_average, C32_average, C33_average, C34_average, C35_average, C36_average, C41_average, C42_average, C43_average, C44_average, C45_average, C46_average, C51_average, C52_average, C53_average, C54_average, C55_average, C56_average,  C61_average, C62_average, C63_average, C64_average, C65_average, C66_average])

        if Lx != 2*d0:
        
            ## Create plot
            fig2 = plt.figure(figsize=(3.5, 2.5))
            plt.subplot(1,1,1)
            ax2 = plt.gca()

            cmap = cm.get_cmap('Blues')
            #cmap = cm.get_cmap('bwr')
            caxes = plt.imshow(C, cmap=cmap)

            Cisotropic= [[r'$\lambda$+2$\mu$', '$\lambda$', '$\lambda$', '$\sim$0', '$\sim$0', '$\sim$0'],['$\lambda$', '$\lambda$+2$\mu$', '$\lambda$', '$\sim$0', '$\sim$0', '$\sim$0'],['$\lambda$', '$\lambda$', '$\lambda$+2$\mu$', '$\sim$0', '$\sim$0', '$\sim$0'], ['$\sim$0', '$\sim$0', '$\sim$0', '$\mu$', '$\sim$0', '$\sim$0'], ['$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\mu$', '$\sim$0'],  ['$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\mu$'] ]

            for i in range(6):
                for j in range(6):
                    c = Cisotropic[i][j]
                    c2 = round(C[i][j],2)
                    if (i==0 and j==0) or (i==1 and j==1) or (i==2 and j==2):               
                        ax2.text(i, j, str(c) + "\n" + str(c2), va='center', ha='center', fontsize = 6, color = 'w')
                    
                    if (i==0 and (j==1 or j==2)) or (i==1 and (j==0 or j==2)) or (i==2 and (j==0 or j==1)) or (i==3 and j==3) or (i==4 and j==4) or (i==5 and j==5):
                        ax2.text(i, j, str(c) + "\n" + str(c2), va='center', ha='center', fontsize = 6, color = 'k')
                    if (i!=j) and ((i==3 or i==4 or i==5) or (j==3 or j==4 or j==5)):
                        ax2.text(i, j, str(c) , va='center', ha='center', fontsize = 6, color = 'k')
        
            plt.xticks([0, 1, 2, 3, 4, 5] , ['C$_{i1}$', 'C$_{i2}$', 'C$_{i3}$', 'C$_{i4}$', 'C$_{i5}$', 'C$_{i6}$'])
            plt.yticks([0, 1, 2, 3, 4, 5] , ['C$_{1j}$', 'C$_{2j}$', 'C$_{3j}$', 'C$_{4j}$', 'C$_{5j}$', 'C$_{6j}$'])
            plt.xticks(fontsize = 8)
            plt.yticks(fontsize = 8)

            # colorbar
            cbar2 = fig2.colorbar(caxes)
            plt.clim(vmax=1)
            cbar2.ax.tick_params(labelsize=8)
            cbar2.set_label(r'GPa', rotation=270, fontsize = 9, labelpad=15)
        
            # title
            plt.title(r'Elastic tensor ' + str(s_title) + r'$d_{0}$', fontsize = 9, pad = 10)
        
            # save plot
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".png")
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".pgf")
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".pdf")
       
        if Lx == 2*d0:
            
            fig2 = plt.figure(figsize=(3.5, 2.5))
            plt.subplot(1,1,1)
            ax2 = plt.gca()

            # plot matrix
        
            #cmap = ListedColormap([ 'lightskyblue', 'steelblue', 'k', 'r'])
            cmap = cm.get_cmap('Blues')
            #cmap = cm.get_cmap('bwr')
            caxes = plt.imshow(C, cmap=cmap)
            
            Cisotropic= [[r'$\lambda$+2$\mu$', '$\lambda$', '$\lambda$', '$\sim$0', '$\sim$0', '$\sim$0'],['$\lambda$', '$\lambda$+2$\mu$', '$\lambda$', '$\sim$0', '$\sim$0', '$\sim$0'],['$\lambda$', '$\lambda$', '$\lambda$+2$\mu$', '$\sim$0', '$\sim$0', '$\sim$0'], ['$\sim$0', '$\sim$0', '$\sim$0', '$\mu$', '$\sim$0', '$\sim$0'], ['$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\mu$', '$\sim$0'],  ['$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\sim$0', '$\mu$'] ]

            for i in range(6):
                for j in range(6):
                    #c = Cisotropic[i][j]
                    c = round(C[i][j],2)
                    if (i==0 and (j==0 or j==1 or j==2)) or (i==1 and (j==0 or j==1 or j==2)) or (i==2 and (j==0 or j==1 or j==2)) or (i==j):               
                        ax2.text(i, j, str(c), va='center', ha='center', fontsize = 6, color = 'w')
                    else:
                        ax2.text(i, j, str(c), va='center', ha='center', fontsize = 6, color = 'k')
            # abscissa
            line1 = ['C$_{i1}$', 'C$_{i2}$', 'C$_{i3}$', 'C$_{i4}$', 'C$_{i5}$', 'C$_{i6}$']
            line2 = ['C$_{1j}$', 'C$_{2j}$', 'C$_{3j}$', 'C$_{4j}$', 'C$_{5j}$', 'C$_{6j}$']
            ax2.set_xticklabels(['']+line1)
            ax2.set_yticklabels(['']+line2)
            plt.xticks(fontsize = 8)
            plt.yticks(fontsize = 8)

            # colorbar
            cbar2 = fig2.colorbar(caxes)
            plt.clim(vmax=1)
            cbar2.ax.tick_params(labelsize=8)
            cbar2.set_label(r'GPa', rotation=270, fontsize = 9, labelpad=15)
        
            # title
            plt.title(r'Elastic tensor $2d_{0}$', fontsize = 9, pad = 10)
        
            # save plot
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".png")
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".pgf")
            fig2.savefig( path_image +  "/Elastic_tensor_" + str(Lx) + ".pdf")
            #pgf_delete_family(path_image2 +  "/Stiffness_matrix_" + str(Lx) + ".pgf", path_image2 +  "/Stiffness_matrix_" + str(Lx) + "_latex.pgf")
    
    #Young Modulus error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()
    
    for l in range(1, len(Size)):
        yerr = [2 * Young_modulus_error_bar[l][0]*10**-9]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Young_modulus_mean[l][0]*10**-9, xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3)
    
    plt.axhline(y=Young_modulus_mean[4][0]*10**-9, color='r', linestyle='--', linewidth = 1)
    ax1.set_ylim(0.732, 1.22)
    ax1.set_yticks([0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20])
    ax1.set_yticks([round(Young_modulus_mean[4][0]*10**-9,3)], minor = True)
    ax1.yaxis.set_ticklabels([round(Young_modulus_mean[4][0]*10**-9,3)], minor = True, fontsize = 8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')
    
    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 10)
    ax1.set_ylabel(r"Young's modulus $E$ (GPa)", fontsize = 9, labelpad = 10)
    
    fig1.savefig( path_image + "/Young_modulus_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Young_modulus_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Young_modulus_error_bar.pdf", bbox_inches="tight")
    plt.show()

    #Bulk Modulus error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    for l in range(1, len(Size)):
        yerr = [2 * Bulk_modulus_error_bar[l][0]*10**-9]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Bulk_modulus_mean[l][0]*10**-9, xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3)

    plt.axhline(y=Bulk_modulus_mean[4][0]*10**-9, color='r', linestyle='--', linewidth = 1)

    ax1.set_ylim(0.495, 0.825)
    ax1.set_yticks([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80])
    ax1.set_yticks([round(Bulk_modulus_mean[4][0]*10**-9,3)], minor = True)
    ax1.yaxis.set_ticklabels([round(Bulk_modulus_mean[4][0]*10**-9,3)], minor = True, fontsize=8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')

    plt.yticks(fontsize = 8) 

    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 10)
    ax1.set_ylabel(r'Bulk modulus $K$ (GPa)', fontsize = 9, labelpad = 10)
    fig1.savefig( path_image + "/Bulk_modulus_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Bulk_modulus_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Bulk_modulus_error_bar.pdf", bbox_inches="tight")
    plt.show()
    
    #Shear Modulus 1 and 2 error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    for l in range(1, len(Size)):
        yerr = [2 * Shear_modulus1_error_bar[l][0]*10**-9]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        yerr2 = [2 * Shear_modulus2_error_bar[l][0]*10**-9]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Shear_modulus1_mean[l][0]*10**-9, xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3, ls='-.')   
        eb2=plt.errorbar(Size_box_mean[l][0]/dmean, Shear_modulus2_mean[l][0]*10**-9, xerr=xerr, fmt='o', color='blue', ecolor='skyblue', elinewidth=1.5, errorevery=2, ls='-.', capsize=3, markersize = 3)
        eb2[-1][0].set_linestyle('--')
        eb3=plt.errorbar(Size_box_mean[l][0]/dmean, Shear_modulus2_mean[l][0]*10**-9, yerr=yerr2, fmt='o', color='blue', ecolor='skyblue', elinewidth=1.5, errorevery=2, ls='-.', capsize=3, markersize = 3)
        eb3[-1][0].set_linestyle('--')
   
        
    plt.axhline(y=Shear_modulus_average_mean[4][0]*10**-9, color='r', linestyle='--', linewidth = 1)
    ax1.set_ylim(0.293, 0.487)
    #ax13.set_yticks([0.24, 0.26, 0.28 ,0.30, 0.32, 0.34, 0.36, 0.38])
    ax1.set_yticks([round(Shear_modulus_average_mean[4][0]*10**-9,3)], minor = True)
    ax1.yaxis.set_ticklabels([round(Shear_modulus_average_mean[4][0]*10**-9,3)], minor = True, fontsize=8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')
    
    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 10)
    ax1.set_ylabel(r'Shear modulus (GPa)', fontsize = 9, labelpad = 10)
    fig1.savefig( path_image + "/Shear_modulus_1_and_2_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Shear_modulus_1_and_2_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Shear_modulus_1_and_2_error_bar.pdf", bbox_inches="tight")
    plt.show()
    
    #Poisson ratio error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    for l in range(1, len(Size)):
        yerr = [2 * Poisson_error_bar[l][0]]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Poisson_mean[l][0], xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3)

    plt.axhline(y=Poisson_mean[4][0], color='r', linestyle='--', linewidth = 1)

    ax1.set_ylim(0, 0.5)
    ax1.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax1.set_yticks([round(Poisson_mean[4][0],3)], minor = True)
    ax1.yaxis.set_ticklabels([round(Poisson_mean[4][0],3)], minor = True, fontsize = 8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')
  
    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 10)
    ax1.set_ylabel(r"Poisson's ratio $\nu$", fontsize = 9, labelpad = 10)
    fig1.savefig( path_image + "/Poisson_ratio_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Poisson_ratio_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Poisson_ratio_error_bar.pdf", bbox_inches="tight")
    plt.show()
    
    #Anisotropy coefficient error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    for l in range(1, len(Size)):
        yerr = [2 * Anisotropy_coefficient_error_bar[l][0]]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Anisotropy_coefficient_mean[l][0], xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3)

    plt.axhline(y=1, color='r', linestyle='--', linewidth = 1)
    ax1.set_ylim(0, 2)
    ax1.set_yticks([0, 0.25, 0.50, 0.75, 1.25, 1.50, 1.75, 2.00])
    ax1.set_yticks([1.00], minor = True)
    ax1.yaxis.set_ticklabels(['1.00'], minor = True, fontsize = 8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')
   
    
    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 9)
    ax1.set_ylabel(r'Zener ratio', fontsize = 9, labelpad = 9)
    fig1.savefig( path_image + "/Anisotropy_coefficient_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Anisotropy_coefficient_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Anisotropy_coefficient_error_bar.pdf", bbox_inches="tight")
    plt.show()
    
    #Density error bar
    fig1 = plt.figure(figsize=(3.5, 2.5))
    plt.subplot(1,1,1)
    ax1 = plt.gca()

    for l in range(1, len(Size)):
        yerr = [2 * Density_error_bar[l][0]]
        xerr = [2 * Length_box_error_bar[l][0]/dmean]
        ax1.errorbar(Size_box_mean[l][0]/dmean, Density_mean[l][0], xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='gray', elinewidth=1.5, capsize=3, markersize = 3)
        
    plt.axhline(y=Density_mean[4][0], color='r', linestyle='--', linewidth = 1)

    ax1.set_ylim(1235, 2058)
    ax1.set_yticks([round(Density_mean[4][0])], minor = True)
    ax1.yaxis.set_ticklabels([round(Density_mean[4][0])], minor = True, fontsize=8)
    ax1.yaxis.set_tick_params(which = 'major', color = 'black', labelcolor = 'black')
    ax1.yaxis.set_tick_params(which = 'minor', color = 'red', labelcolor = 'red')   
        
        
    plt.xticks(np.arange(0, ceil(Size_box_mean[len(Size)-1][0]/dmean)+1, 2.0))  
    plt.xticks(fontsize = 8)
    plt.yticks(fontsize = 8)
    ax1.set_xlabel(r'$L_{x}/d_{0}$', fontsize = 9, labelpad = 10)
    ax1.set_ylabel(r'density (Kg/m$^3$)', fontsize = 9, labelpad = 10)
    fig1.savefig( path_image + "/Density_error_bar.png", bbox_inches="tight")
    fig1.savefig( path_image + "/Density_error_bar.pgf", bbox_inches="tight")
    fig1.savefig( path_image + "/Density_error_bar.pdf", bbox_inches="tight")
    plt.show()
    
 
