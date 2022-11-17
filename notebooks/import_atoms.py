from Atom_class import Atom

def import_atoms(filename):
    atoms = []

    file = open(filename, "r")
    current_read = "header"

    for line in file:
        if line == "\n":
            continue

        if current_read == "header":
            if "Atoms" in line:
                current_read = "atoms"
                continue
            
            elif "xlo" in line:
                split = line.split(" ")
                xlo = float(split[0])
                xhi = float(split[1])
            
            elif "ylo" in line:
                split = line.split(" ")
                ylo = float(split[0])
                yhi = float(split[1])
            
            elif "zlo" in line:
                split = line.split(" ")
                zlo = float(split[0])
                zhi = float(split[1])
            
        elif current_read == "atoms":
            if "Velocities" in line:
                current_read = "velocities"
                continue

            atoms.append(Atom(*line.split(" ")[:7]))

        else: # current_read == "velocities"
            pass

    file.close()

    boundaries = [xlo, xhi, ylo, yhi, zlo, zhi]
    return atoms, boundaries

