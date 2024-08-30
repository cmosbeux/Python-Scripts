#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:36:03 2022

@author: cmosbeux

Description: Functions for data plot

"""
from MyColobars import alpha_cmap
import numpy as np
import matplotlib.cm as cm
from format_reading import netcdf
import copy
import matplotlib as mpl
import h5py


#%%
#Get some observation plot
def plot_velocity_map(ax):
    
    #load data
    x, y, vx = netcdf.readgrid('/Users/cmosbeux/Documents/Data/antarctica_velocity_updated_v2.nc', 'vx')
    x, y, vy = netcdf.readgrid('/Users/cmosbeux/Documents/Data/antarctica_velocity_updated_v2.nc', 'vy')
    v = (vx**2+vy**2)**0.5
    
    #make cmap and plot figure
    cmap = cm.get_cmap('Blues', 21)  
    cb = ax.imshow(v.data, cmap = cmap, extent = [x[0], x[-1], y[-1], y[0]], norm = mpl.colors.LogNorm(vmin=1, vmax=1000), interpolation = 'bilinear', alpha=0.75)
    
    #colorbar
    # cbar = ax.cax.colorbar(cb)
    # cbar = ax.cbar_axes[0].colorbar(cb)
    # cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015]) 
    # cbar = fig.colorbar(cb, cax=cbaxes, ticks=[1,100,1000],orientation='horizontal')
    # cbar.set_label('ice surface speed', fontsize=16)
    # cbar.ax.set_xticklabels(['1', '100', r'   1000 m a$^{-1}$'], fontsize=14)
    
    
def plot_melt_obs(ax):
    #load data
    filename  = '/Users/cmosbeux/Documents/Data/ANT_iceshelf_melt_rates_CS2_2010-2018_v0_5km.h5'
    f = h5py.File(filename, "r")
    # List all groups
    #print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    x = f['x']
    y = f['y']
    data = f[a_group_key]
    
    #make cmap and plot figure
    cmap = cm.get_cmap('hot_r', 21)  
    cb = ax.imshow(data[()][:,:].T, cmap = cmap, extent = [x[0][0], x[0][-1], y[0][-1], y[0][0]], norm = mpl.colors.Normalize(vmin=0, vmax=20))
 
    
def plot_dhdt_obs(ax, shelf=True):
    #load data
    filename  = '/Users/cmosbeux/Documents/Data/ICESat/AIS_mass_change.h5'
    
    f = h5py.File(filename, "r")
    # List all groups
    #print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    # Get the data
    x = f['x']
    y = f['y']
    data = f['dHdt']
    #make cmap and plot figure
    cmap = cm.get_cmap('seismic', 21)  
    dx = 7e4
    dy = 2.5e4
    data = data[()][:,::-1].T
    #cb = ax.imshow(, cmap = cmap,  extent = [-2.595e6+dx, 2.58e6+dx, -2.16e6+dy, 2.16e6+dy], norm = mpl.colors.Normalize(vmin=-0.1, vmax=0.1))  
    vmin, vmax = -10,10
    data_neg = copy.copy(data)
    data_neg[data_neg>0]=np.nan
    cmap = alpha_cmap(cm.seismic, type = 'max', width=128)
    cax2 = ax.imshow(data_neg, cmap = cmap,  extent = [-2.595e6+dx, 2.58e6+dx, -2.16e6+dy, 2.16e6+dy], norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    data_pos = copy.copy(data)
    data_pos[data_pos<0]=np.nan 
    vmin, vmax = -1,  1
    cmap = alpha_cmap(cm.seismic, type = 'min', width=128)
    cax1= ax.imshow(data_pos, cmap = cmap,  extent = [-2.595e6+dx, 2.58e6+dx, -2.16e6+dy, 2.16e6+dy], norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    
    
    #Shelf
    if shelf:
        filename = '/Users/cmosbeux/Documents/Data/ICESat/ICE1_ICE2_AnIS_dHdt_2003_2018_R209_05KM_FLOAT_MASS_F2.h5'
        f = h5py.File(filename, "r")
        # List all groups
        #print("Keys: %s" % f.keys())
        a_group_key = list(f.keys())[0]
        # Get the data
        x = f['x']
        y = f['y']
        data = f['dhdt']
        #make cmap and plot figure
        cmap = cm.get_cmap('seismic', 21)  
        dx = 7e4
        dy = 2.5e4
        data = data[()]*9.26
        #cb = ax.imshow(, cmap = cmap,  extent = [-2.595e6+dx, 2.58e6+dx, -2.16e6+dy, 2.16e6+dy], norm = mpl.colors.Normalize(vmin=-0.1, vmax=0.1))  
        vmin, vmax = -10, 10
        data_neg = copy.copy(data)
        data_neg[data_neg>0]=np.nan
        cmap = alpha_cmap(cm.seismic, type = 'max', width=128)     
        extent = [np.min(f['x']), np.max(f['x']), np.min(f['y']), np.max(f['y'])]
        cax2 = ax.imshow(data_neg, cmap = cmap, extent = extent , norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax))
        data_pos = copy.copy(data)
        data_pos[data_pos<0]=np.nan
        vmin, vmax = -1,  1
        cmap = alpha_cmap(cm.seismic, type = 'min', width=128)
        cax1= ax.imshow(data_pos, cmap = cmap,  extent = extent, norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    
