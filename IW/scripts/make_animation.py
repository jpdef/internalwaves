import feather
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import functools
import cmocean
import glob

def plot_depth_slice(df,ax,index,pcol):
    x = df['x'].unique()
    y = df['y'].unique()
    z = df['z'].unique()[index]
    xx,yy =np.meshgrid(x,y)
    zeta = np.array(df.loc[df['z']==z][pcol])
    zeta = np.reshape(zeta,(len(x),len(y)))
    p = ax.contourf(xx/1000,yy/1000,zeta,20,cmap=cmocean.cm.thermal)
    ax.set_xlabel("Latitude (km)")
    ax.set_ylabel("Longitude (km)")
    return(p)
    
def plot_lon_slice(df,ax,index,pcol):
    x = df['x'].unique()
    z = df['z'].unique()
    y = df['y'].unique()[index]
    xx,zz =np.meshgrid(x,z)
    zeta = np.array(df.loc[df['y']==y][pcol])
    zeta = np.reshape(zeta,(len(x),len(z)))
    p = ax.contourf(xx/1000,zz/1000,zeta,20,cmap=cmocean.cm.thermal)
    ax.set_xlabel("Latitude (km)")
    ax.set_ylabel("Depth (km)")
    return p

def plot_lat_slice(df,ax,index,pcol):
    y = df['y'].unique()
    z = df['z'].unique()
    x = df['x'].unique()[index]
    yy,zz =np.meshgrid(y,z)
    zeta = np.array(df.loc[df['x']==x][pcol])
    zeta = np.reshape(zeta,(len(y),len(z)))
    p = ax.contourf(yy/1000,zz/1000,zeta,20,cmap=cmocean.cm.thermal)
    ax.set_xlabel("Longitude (km)")
    ax.set_ylabel("Depth (km)")
    return p

 
def make_animation(files,plot_func):
    fig, ax = plt.subplots()
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    set_animation_attributes(fig,ax) 
    
    #First Frame
    df = feather.read_dataframe(files[0])
    p = plot_func(df,ax)
    
    #Init and update function
    init = functools.partial(init_animation,fig,ax,p,cbar_ax)
    update = functools.partial(update_animation,fig,ax,cbar_ax,plot_func)
   
    #Compile Animation
    print("Creating Animation it may take awhile")
    ani = FuncAnimation(fig, update, frames=files[1:],init_func=init, blit=False)
    path = 'animation.mp4'
    ani.save(path)
    return ani


def set_animation_attributes(fig,ax):
    ax.invert_yaxis()
    ax.set_xlabel('Range Km')
    ax.set_ylabel('Depth Km')

    
def init_animation(fig,ax,p,cbar_ax):
    fig.subplots_adjust(right=0.8)
    cbar = fig.colorbar(p, cax=cbar_ax)
    cbar.set_label('Sound Speed Anom. (m/s)')
    fig.set_dpi(300)
    return p


def update_animation(fig,ax,cbar_ax,plot_func,frame):
    
    df = feather.read_dataframe(frame)
    ax.set_title('%.2f hours' % (df['t'][0]/3600.0))
    
    p = plot_func(df,ax)
    cbar_ax.cla()
    fig.colorbar(p, cax=cbar_ax)
    return p

