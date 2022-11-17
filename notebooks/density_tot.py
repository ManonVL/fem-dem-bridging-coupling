import math
import numpy as np
from Atom_class import Atom
from import_atoms import import_atoms


def density_tot(input_file):
    
    atoms, boundaries = import_atoms(input_file)
    
    Lx = abs(boundaries[1] - boundaries[0])
    Ly = abs(boundaries[3] - boundaries[2])
    Lz = abs(boundaries[5] - boundaries[4])
    
    mp = 0
    
    for atom in atoms:
        mp = mp + 4/3*math.pi*(float(atom.diameter)/2)**3 * atom.density


    V = Lx * Ly * Lz
    density_tot = mp/V

    return(density_tot)



