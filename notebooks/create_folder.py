import os

def create_folder(path):
    try:
        os.mkdir(path)
    except OSError:
        print ("")
    else:
        print ("Successfully created the directory %s" % path)
