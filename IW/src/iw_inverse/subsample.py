import glob
import os

def get_file_list(path,fpat='run*'):
    dl = glob.glob(os.path.join(path,fpat))
    #print(dl)
    #dl = sorted(dl)
    return dl

#def read_data_dir(path):
#    for f in get_file_list(path):
#        pass

