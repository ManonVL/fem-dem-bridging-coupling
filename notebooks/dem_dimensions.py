def dem_dimensions(path, file_name):
    
    open1 =  path + "/" + file_name
    f1 = open(open1, "r")
    Adata = []
    Adata = [ line.split() for line in f1]
    DEMLX = abs(float(Adata[5][1])-float(Adata[5][0]))
    DEMLY = abs(float(Adata[6][1])-float(Adata[6][0]))
    DEMLZ = abs(float(Adata[7][1])-float(Adata[7][0]))

    return DEMLX, DEMLY, DEMLZ
