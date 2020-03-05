import os
import feather
import iw_misc
from tqdm import tqdm

#Maps a new columns onto a set of dataframes
def map_scalars(path,out_col,fn):
   files = os.listdir(path)
   for f in tqdm(files,ascii=True,
                 total=len(files),
                 leave=True,desc="Mapping Scalars"):
        fpath = os.path.join(path,f)
        df = feather.read_dataframe(fpath)
        if 'd' in df.columns:
           df[out_col] = fn(df)
        feather.write_dataframe(df,fpath)


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


