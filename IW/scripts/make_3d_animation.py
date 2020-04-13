import feather
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
import numpy as np
import functools
import cmocean
import glob
import os

def plot_3d_cube(df,ax,index,pcol):
    cmap = mpl.cm.RdBu
    #cmap = cmocean.cm.thermal
    
    dmax = max(df[pcol].unique())
    dmin = min(df[pcol].unique())
    davg = np.mean(df[pcol].unique())
     
    
    norm = mpl.colors.DivergingNorm(vmin=dmin,vcenter=davg, vmax=dmax)

    x = df['x'].unique()
    y = df['y'].unique()
    z = df['z'].unique()
    
    xx,yy =np.meshgrid(x,y)
    xx,zz =np.meshgrid(x,z)
    yt,zt =np.meshgrid(y,z)
    
    #TOP FACE
    zindex = z[index]
    zeta = np.array(df.loc[df['z']==zindex][pcol])
    zeta = np.reshape(zeta,(len(x),len(y)))
    p = ax.contourf(xx/1000,yy/1000,zeta,20,zdir='z', offset=0,
                                    cmap=cmap,norm=norm)
    
    #LON SLICE
    yindex = y[index]
    zeta = np.array(df.loc[df['y']==yindex][pcol])
    zeta = np.reshape(zeta,(len(x),len(z)))
    p = ax.contourf( xx/1000,zeta,zz/1000,20,zdir='y', offset=0,
                                    cmap=cmap,norm=norm)
    
    #LAT SLICE
    xindex = x[index]
    zeta = np.array(df.loc[df['x']==xindex][pcol])
    zeta = np.reshape(zeta,(len(y),len(z)))
    p = ax.contourf(zeta,yt/1000,zt/1000,20, zdir='x', offset=0,
                                    cmap=cmap,norm=norm)
    
    return(p)
    
from mpl_toolkits.mplot3d import Axes3D 
def make_animation(path,plot_func,findex=0):

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)
    set_animation_attributes(fig,ax) 
    
    files = get_file_list(path)
    files = files[:findex] if findex > 0 else files 
     
    #First Frame
    df = feather.read_dataframe(files[0])
    scale_axis(df,ax) 
    p = plot_func(df,ax)
    
    #Init and update function
    init = functools.partial(init_animation,fig,ax,p,cbar_ax)
    update = functools.partial(update_animation,fig,ax,text,cbar_ax,plot_func)
   
    #Compile Animation
    print("Creating Animation it may take awhile")
    ani = FuncAnimation(fig, update, frames=files[1:],init_func=init, blit=False)
    path = 'animation.mp4'
    ani.save(path)
    return ani


def scale_axis(df,ax):
    xmax = max(df['x'].unique()/1000)
    xmin = min(df['x'].unique()/1000)
    
    ymax = max(df['y'].unique()/1000)
    ymin = min(df['y'].unique()/1000)
    
    zmax = max(df['z'].unique()/1000)
    zmin = min(df['z'].unique()/1000)

    ax.set_xlim3d(xmax,xmin)
    ax.set_ylim3d(ymin,ymax)
    ax.set_zlim3d(zmax,zmin)


def set_animation_attributes(fig,ax):
    ax.set_xlabel('Range Km')
    ax.set_zlabel('Depth Km')


    
def init_animation(fig,ax,p,cbar_ax):
    fig.subplots_adjust(right=0.8)
    cbar = fig.colorbar(p,cax=cbar_ax)
    cbar.set_label('Sound Speed Anom. (m/s)')
    fig.set_dpi(300)
    return p


def update_animation(fig,ax,text,cbar_ax,plot_func,frame):
    
    df = feather.read_dataframe(frame)
    text.set_text( '%.2f hours' % (df['t'][0]/3600.0) )
    p = plot_func(df,ax)
    #cbar_ax.cla()
    #fig.colorbar(p, cax=cbar_ax)
    return p

def get_file_list(path,fpat="run*"):
    files = glob.glob(os.path.join(path,fpat))
    return sorted(files)
