import os
import feather
import iw_misc
from tqdm import tqdm

#Maps a new columns onto a set of dataframes
def map_scalars(path,out_col,fn):
   files = os.listdir(path)
   for f in tqdm(files,ascii=True,
                 total=len(files),
                 leave=True,desc="Mapping Scalars -> " + out_col):
        fpath = os.path.join(path,f)
        df = feather.read_dataframe(fpath)
        if 'd' in df.columns:
           df[out_col] = fn(df)
        feather.write_dataframe(df,fpath)

def map_anom(df,mean,grad,in1='z',in2='d'):
    z = df[in1]
    d = df[in2]
    return mean(z) + grad(z)*d

def map_svel(df,func):
    return func(df['S'],df['T'],df['z'])

def map_svel_anom(df,func):
    return df['CC'] - func(df['z'])


#Exampling mapping function for sound speed
def map_sound_speed(df,ssp=None,ssp_grad=None):
    z = df['z']
    d = df['d']
    if ssp and ssp_grad:
        return ssp(z) + ssp_grad(z)*d
    else:
        return iw_misc.sound_prof_munk(z) + iw_misc.sound_grad_munk(z)*d


def map_sound_speed_anom(df,ssp_grad=None):
    z = df['z']
    d = df['d']
    
    if ssp_grad:
        return ssp_grad(z)*d
    else:
        return iw_misc.sound_grad_munk(z)*d


