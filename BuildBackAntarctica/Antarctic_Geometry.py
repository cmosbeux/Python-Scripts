#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 10:04:19 2023

@author: cmosbeux

Description: This script reconstructs the ice sheet as it was back in time based on current geometry and dhdt trend. The grounding line position and flotation criterion 
             are based on hydrostatic equilibrium hypothesis. 
             The reconstrution for a given year is based on the hypothesis that bedmachine geometry is representative  of the year 2020.

Requierment: 
        - topography file (containing ice thickness and bedrock geometry): e.g.: BedMachineAntarctica-v03.nc
        - dhdt files for the shelves and the ice sheet, I uses dhdt observations from IceSat1-2 missions over 2003-2019

Potential Requierment:
        - Path to data are hard coded for my system, you should define the filename in the different functions
        - if you want to make some plots, you will need my Antarctica_Background Module and add it to your paths (available on my github)
        - I handle netcdf with a format_reading Module of my own too (also available on my github)
        - You might want to adjust water and ice density depending on the densities you will use in your modelling (ri and rw).
"""


from Antarctica_Background import Plot_Antarctica, basin, scale,plot_GL_MEASURES
from format_reading import netcdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import cm
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import h5py
from FunctionsDataPlot import plot_dhdt_obs
import copy
from scipy.ndimage import gaussian_filter
import sys
from scipy import interpolate
import os

#-------------------------------------------------------------------------
# USER KEYWORDS
#-------------------------------------------------------------------------
plot = True
ri = 0.917
rw = 1.028
z_sl = 0
targetyear = 1850
resolution = 1 #km
power_rate = 1 # for long period, we can use dhdt rates that decrease with t
#-------------------------------------------------------------------------

#%%Read dHdT from IceSat1-2 (over the shelves)
def Build_IceSat2_ContinuousMap(x_final=None, y_final=None, masked=False):
    print("Building IceSat dhdt...")
    #Ice Sheet
    filename  = '/Users/cmosbeux/Documents/Data/ICESat/AIS_mass_change.h5'
    
    f = h5py.File(filename, "r")
    
    a_group_key = list(f.keys())[0]
    # Get the data
    x_ground = f['x'][()].T
    y_ground = f['y'][()].T[::-1,:]
    
    data_ground = f['dHdt'][()].T[::-1,:]
    
    #Ice Shelf
    filename = '/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/ICESat2/ICE1_ICE2_AnIS_dHdt_2003_2018_R209_05KM_FLOAT_MASS_F2.h5'
    f = h5py.File(filename, "r")
    a_group_key = list(f.keys())[0]
    # Get the data
    x_shelf = f['x'][()]
    y_shelf = f['y'][()]
    data_shelf = f['dhdt'][()]
    data_shelf = data_shelf*9.26
    
    #Interpolation of the shelf data on the ground data
    print("\tInterpolating over same grid...")
    lim = 1
    x_share, y_share = x_ground[lim:-lim,lim:-lim], y_ground[lim:-lim,lim:-lim] 
    # print('xground = %d' % x_share.min())
    points_share = list(zip(x_share.flatten(), y_share.flatten()))
    
    #make the interpolator
    interpolator = RegularGridInterpolator((x_shelf[0,:], y_shelf[::-1,0]), data_shelf[::-1,:].T)
    #interpolation
    data_shelf_share = interpolator(points_share)
    data_shelf = data_shelf_share.reshape(x_share.shape)
    
    #mask 
    mask = np.logical_and(np.isnan(data_shelf), np.isnan(data_ground[lim:-lim,lim:-lim]))
    n, m = mask.shape
    
    mask2 = copy.copy(mask)
    for i in np.arange(n):
        for j in np.arange(m):
            if mask[i,j]:
                if False in mask[i-1:i+2,j-1:j+2]:
                    mask2[i,j] = False
                    
    data = np.nan_to_num(data_shelf) + np.nan_to_num(data_ground[lim:-lim,lim:-lim])
    data[data==0.0]=np.nan
    
    #mask NaNs
    data = np.ma.masked_invalid(data)
    
    #get only the valid values
    x1 = x_share[~data.mask]
    y1 = y_share[~data.mask]
    new_data = data[~data.mask]
    
    new_data = interpolate.griddata((x1, y1), new_data.ravel(), (x_share, y_share),method='linear')
    data = np.ma.MaskedArray(new_data, mask2)
    
    #reshape to a given output grid
    if (x_final is not None) and (y_final is not None):
        print("\tRegridding")
        print("\t - origin boundaries : x = [%d, %d], y = [%d, %d]" % (x_share.min(), x_share.max(), y_share.min(), y_share.max()))
        print("\t - target boundaries : x = [%d, %d], y = [%d, %d]" % (x_final.min(), x_final.max(), y_final.min(), y_final.max()))
        
        y_offset = 5e4 #I have to apply offset but I don't know why
        
        interpolator2 = RegularGridInterpolator((x_share[0,:], y_share[::-1,0] - y_offset), data.T, bounds_error=False)
        xx,yy = np.meshgrid(x_final,y_final)
        points_final = list(zip(xx.flatten(), yy.flatten()))
        
        data = interpolator2(points_final)
        data = data.reshape((len(x_final),len(y_final)))
        x, y = x_final, y_final
        
    else: 
        x, y = x_share, y_share
    print("\tDone.")
    return x, y, data


def read_bedmachine(resolution=5):
    print('Reading bedmachine...')
    print('\t -using a %0.1f km resolution' % (resolution*0.5))
    bedmachine_filename = '/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/BedMachineAntarctica-v03.nc'
    
    x, y, bed = netcdf.readgrid(bedmachine_filename, 'bed')
    x, y, surf = netcdf.readgrid(bedmachine_filename, 'surface')
    x, y, thickness = netcdf.readgrid(bedmachine_filename, 'thickness')
    x, y, mask = netcdf.readgrid(bedmachine_filename, 'mask')
    
    x_inf, x_sup = +0, -1
    y_inf, y_sup = +0, -1
    x = x[x_inf:x_sup:resolution]
    y = y[y_inf:y_sup:resolution]
    bed = bed[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
    surf = surf[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
    thickness = thickness[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
    mask=mask[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
    
    mask2 = np.zeros_like(mask)
    mask2[mask==0] = True
    mask2[mask>0.5] = False
    
    bed = np.ma.masked_array(bed,mask2)
    surf = np.ma.masked_array(surf,mask2)
    thickness = np.ma.masked_array(thickness,mask2)
    
    version = bedmachine_filename.split('/')[-1][-6:-3]
    
    return x,y,bed,surf,thickness, version


def grounding(thickness, bed, ri=0.917, rw = 1.028, z_sl=0):
    cond = thickness*ri/rw + bed
    gmask = np.ones_like(cond)
    gmask[cond<0] = -1  #floating
    gmask[cond>0] = 1 #grounded
    return gmask


def new_floating_surf(thickness, bed, gmask, ri=0.917, rw = 1.028, zs_sl=0):
    surface_ground = bed + thickness
    surface_float = thickness*(1-ri/rw)
    grounded_mask = 1+((gmask-1)/2) 
    floating_mask = (-((gmask-1)/2))
    surface = surface_ground*grounded_mask + surface_float*floating_mask
    return surface


def get_GL(x, y, gmask):
    print('Extracting grounding line...')
    x_gl,y_gl = [],[]
    n, m = len(x), len(y)
    for i in range(n):
        for j in range(m):
            if gmask[i,j] == 1 and np.mean(gmask[i-1:i+2,j-1:j+2])!=1:
                x_gl.append(x[j])
                y_gl.append(y[i])
        old_i = 0
        if int((i/n+1e-2)*100) != old_i:
            sys.stdout.write('\r')
            sys.stdout.flush()
            sys.stdout.write(" [%d%%]" % ((i/n+1e-2)*100))
            old_i = int((i/n+1e-2)*100)
            sys.stdout.flush()
    print('\n')
    return x_gl, y_gl


def create_directory(directory_path):
    """
    Create a directory if it does not exist.

    Parameters:
    - directory_path (str): The path of the directory to be created.
    """
    try:
        # Check if the directory already exists
        if not os.path.exists(directory_path):
            # Create the directory
            os.makedirs(directory_path)
            print(f"Directory '{directory_path}' created successfully.")
        else:
            print(f"Directory '{directory_path}' already exists.")
    except Exception as e:
        print(f"Error creating directory '{directory_path}': {str(e)}")

        
def compute_rate(year0, year1, power_rate):
    annual_rate = []  
    years = np.arange(year0, year1,1)
    for i,year in enumerate(years[:-1]):
        annual_rate.append(((years[i]-years[0])/(years[-1]-years[0]))*power_rate)
    rate = sum(annual_rate)/(years[-1]-years[0])
    return rate

#%% Loading files
x, y, bed, surf, thickness, version = read_bedmachine(resolution=resolution*2)
x_IS, y_IS, dhdt_IS = Build_IceSat2_ContinuousMap(x_final=x, y_final=y)

#make directory

create_directory(f'{targetyear}_State')


#%% build dhdt trend and corresponding dh
#the trend is for the period 2003-2019, we assume it is similar over our period of reconstruction
# extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)] 

extent_IS = [np.min(x), np.max(x), np.min(y), np.max(y)]

dhdt_IS = dhdt_IS[::-1,:]
dhdt_IS = np.ma.masked_array(dhdt_IS.data, surf.mask)


rate = compute_rate(targetyear, 2022, power_rate)
dh_IS_targetyear = dhdt_IS*(2022-targetyear)*rate 
dh_IS_2020 = dhdt_IS*(2022-2020)


base = surf - thickness

gmask = np.ones_like(base)
gmask[base-bed>0] = -1

# target year
thickness_targetyear = thickness - np.nan_to_num(dh_IS_targetyear)
gmask_targetyear = grounding(thickness_targetyear, bed)
surf_targetyear = new_floating_surf(thickness_targetyear, bed, gmask_targetyear)
bed_targetyear = surf_targetyear - thickness_targetyear
bedrock_targetyear = bed

#2020
thickness_2020 = thickness - np.nan_to_num(dh_IS_2020)
gmask_2020 = grounding(thickness_2020, bed)
surf_2020 = new_floating_surf(thickness_2020, bed, gmask_2020)
bed_2020 = surf_2020 - thickness_2020
bedrock_2020 = bed

#%%
thickness_targetyear[thickness_targetyear<10]=10.
thickness_targetyear[thickness_targetyear<10]=10.

#%% Get grounding lines for target year and 2020
GL = True
if GL:
    x_gl_2020, y_gl_2020 = get_GL(x,y,gmask_2020) 
    x_gl_targetyear, y_gl_targetyear = get_GL(x,y,gmask_targetyear) 
    np.savetxt('2020_State/GL_2020.txt', [x_gl_2020, y_gl_2020])
    np.savetxt(f'{targetyear}_State/GL_{targetyear}.txt', [x_gl_targetyear, y_gl_targetyear])
else:
    x_gl_2020, y_gl_2020 = np.loadtxt('2020_State/GL_2020.txt')
    x_gl_targetyear, y_gl_targetyear = np.loadtxt(f'{targetyear}_State/GL_{targetyear}.txt')
#%%
from Antarctica_Background import Plot_Antarctica, basin, scale, plot_front, basemap_LIMA_AMU


if plot:
    extent = np.asarray(basin.PanAntarctic())
    fig, ax = Plot_Antarctica(nrows=1, ncols=1, basemap = 'light', GL = False, icefront = False, continental_shelf=0.0, extent=extent, figsize = (30,25))
    #basemap_LIMA_AMU(ax[0])
    cb = ax[0].imshow(dhdt_IS*rate, cmap = 'seismic',extent=extent_IS, vmin = -10, vmax = 10, alpha = 0.9, zorder = 1e1)
    ax[0].scatter(x_gl_2020, y_gl_2020, c='black', s=0.1, linewidths=0, zorder= 1e5)
    ax[0].scatter(x_gl_targetyear, y_gl_targetyear, c='orange', s=0.1, linewidths=0, zorder = 1e5)
    
    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.35, 0.12, 0.013]) 
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[-10,0,10],orientation='horizontal')
    cbar.set_label('ice thickness rate of change', fontsize=16)
    cbar.ax.set_xticklabels(['-10', '0', r'   10 m a$^{-1}$'], fontsize=13)   
    ax[0].axis('off')

    plot_GL_MEASURES(ax[0], '1996', color = 'green', lw = 0.3, zorder = 1e5)
    
    plt.savefig(f'Figures/Antarctica_dHdt_GL{targetyear}.pdf', bbox_inches='tight', pad_inches=0.1, transparent = True)
    plt.close('all')
#%%
#targetyear State
dic = {'x':x, 'y':y, 'surface':surf_targetyear, 'bed':bedrock_targetyear, 'thickness':thickness_targetyear}
netcdf.write(f'{targetyear}_State/BedMachine_IS_{targetyear}_{version}.nc' , dic)

#2020 State
dic = {'x':x, 'y':y, 'surface':surf_2020, 'bed':bedrock_2020, 'thickness':thickness_2020}
netcdf.write('2020_State/BedMachine_IS_2020_{version}.nc', dic)