#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:07:19 2020

@author: cmosbeux
"""

import os
from format_reading import netcdf
import numpy as np

def map_choice(name, variable):
    
    """choice the map and the associated variable (e.g. bedrock, surface, base, thickness, cavities)."""
    
    map_names = ['bedmap2', 'bedmachine']
    
    if name == 'bedmap2':
        path = os.path.expanduser('~/Documents/Developper/Data/bedmap2.nc')
        if variable == 'base': variable = 'bottom'
        #bedmap2 (x and y are reversed in bedmap2 vs bedmachine)
        y, x, z = netcdf.readgrid(path, variable )

    elif name == 'ris_melt':
        path = os.path.expanduser('~/Documents/Developper/Data/RIS_Melt_Monthly.nc')
        if variable == 'meltrate':
            variable = 'BMB_mean'
        x, y, z = netcdf.readgrid(path, variable )                        
    
    elif name == 'ris_bedmap2':
        path = os.path.expanduser('~/Documents/Developper/Data/ros_bedmap2.nc')
        if variable == 'base': variable = 'bottom'
        #bedmap2 (x and y are reversed in bedmap2 vs bedmachine)
        y, x, z = netcdf.readgrid(path, variable )
        
    elif name == 'bedmachine':
        path = os.path.expanduser('~/Documents/Developper/Data/BedMachineAntarctica_2019-11-05_v01.nc')
        if variable == 'bedrock': variable = 'bed'
        if variable == 'base':
            x, y, surface = netcdf.readgrid(path, 'surface')
            x, y, thickness = netcdf.readgrid(path, 'thickness')
            z = surface - thickness
        else:
            y, x, z = netcdf.readgrid(path, variable )
            
  
    else:
        print('no such map in database. List of map:')
        for i in map_names:
            print('\t- %s' % i)
            
    return x,y,z
#%%
    
def map_delimitation(limit, x,y, Data, eps =1e-3):
    """Delimit a subset of a 2D array with limit = [xinf,xsup, yinf, ysup]"""
    X,Y = np.meshgrid(x,y)
    Xspan = np.where((X < limit[1]+eps) & (X > limit[0]-eps))[1][[0, -1]]
    Yspan = np.where((Y < limit[3]+eps) & (Y > limit[2]-eps))[0][[0, -1]]
    print(Xspan, Yspan)
    # Create a selection
    sel = [slice(Xspan[0], Xspan[1] + 1), slice(Yspan[0], Yspan[1] + 1)]
    sel = tuple(sel)
    
    # Extract
    newX = X[sel]  # == Lons[Xspan[0]:Xspan[1]+1, Yspan[0]:Yspan[1]+1]
    newY = Y[sel]
    newData = Data[sel]
    
    return newX, newY, newData



# limit = [40, 60, 20, 50]
# print limit
# x = range(0,100)
# y = range(0,100)          
# Data = np.random.randn(100,100)            
# xdel, ydel, zdel = map_delimitation(limit,x,y,Data)

            
            