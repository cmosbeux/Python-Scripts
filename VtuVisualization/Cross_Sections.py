#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:09:47 2022

@author: cmosbeux
"""

from FunctionsVTU import vtu2section
from Antarctica_Background import basin, Plot_Antarctica, scale, colorbar, Plot_Antarctica2
import numpy as np
import matplotlib.pyplot as plt
import ParticleTracking

'''Give the starting position of the flowline'''
#PIG
section = (-1.57e6, -1.5e5)
#Crosson 1
#section = (-1470888, -620738)
#section = (-1.45e6, -5.82e5)

#Crosson 2
# x0,y0 = -1.54e6, -5.65e5
# x1,y1 = -1.4587e6, -6.635e5   
# section = [x0,y0,x,y1]
# section = [-1546902, -546342,-1470888, -620738]

map_extent = np.asarray(basin.Amundsen())
#%%

x,y,zs,zb,bed = [],[],[],[],[]
dist = []
time = []

periods = ('0-10', '20-30', '20-30', '50-60','90-100')
step = (1,1,11,1,18)

fixed_flowline = False #compute the flowline at first time step
for s,p in zip(periods,step): 
    filename = '../Output@dahu/vtuoutputs/PICO_GL/amundsen_s6_pico_tp3_%sy_approx__t%0.4d.pvtu' % (s,p)
    x_i,y_i,zs_i, zb_i, bed_i = vtu2section(filename, map_extent, section, length=320, dx=1e3, fixed_flowline=fixed_flowline)
    x_i,y_i,zs_i,zb_i,bed_i = np.asarray(x_i), np.asarray(y_i), np.asarray(zs_i), np.asarray(zb_i), np.asarray(bed_i)
    fixed_flowline = True
    zs_i[zs_i<zb_i]=np.nan
    zs_i[zs_i<0]=-0.1*zb_i[zs_i<0]
    zb_i[np.isnan(zs_i)]=np.nan
    
    zb_i[zb_i<bed_i]=bed_i[zb_i<bed_i]
    
    #check first nan
    if np.isnan(zs_i[0]):
        case = 1
        for i in np.arange(len(zs_i)-1,-1,-1):
            if np.isnan(zs_i[i]):
                break
            else:
                continue
    else:
        case = 2
        for i in np.arange(0,len(zs_i),1):
            if np.isnan(zs_i[i]):
                break
            else:
                continue
        
    #distance to front
    dist_i = ((x_i-x_i[0])**2+(y_i-y_i[0])**2)**0.5/1e3

    #to get a vertical front in the plot (only visual)
    if i!=0 and case == 1:
        zs_i[i]=zb_i[i+1]
        zs_i[:i-2]=np.nan
        zb_i[:i-2]=np.nan
        dist_i[i+1]=dist_i[i]
    elif (i!=len(dist_i)-1) and case == 2:
        zs_i[i]=zb_i[i-1]
        zs_i[i+2:]=np.nan
        zb_i[i+2:]=np.nan
        dist_i[i-1]=dist_i[i]
    elif (i==len(dist_i)-1) and case == 2:
        zs_i[i]=zb_i[i-1]
        dist_i[i-1]=dist_i[i]
        
        
    #build the lists
    dist.append(dist_i)
    x.append(x_i),y.append(y_i),zs.append(zs_i), zb.append(zb_i), bed.append(bed_i)
    time.append(s.split('-')[0])

#%% Plot
#-------

#trick to get the front to dist = 0
for i in range(0,len(dist)):
    dist[i] = dist[i].max()-dist[i]

#bed
fig1, ax = plt.subplots(1,1, figsize=(20,10))
plt.plot(dist[0], bed[0], 'tan', zorder = 1000 )

#surface
ax.plot(dist[0], zs[0], 'k--', label='initial' )
ax.plot(dist[1], zs[1], 'red', label='20 years' )
ax.plot(dist[2], zs[2], 'turquoise', label='30 years' )
ax.plot(dist[3], zs[3], 'lightseagreen', label='50 years' )
ax.plot(dist[4], zs[4], 'coral', label='100 years' )

#base
ax.plot(dist[0], zb[0], 'k--', zorder=-10 )
ax.plot(dist[1], zb[1], 'red' , zorder=-10)
ax.plot(dist[2], zb[2], 'turquoise', ls='-', zorder = -20)
ax.plot(dist[3], zb[3], 'lightseagreen', zorder=-10)
ax.plot(dist[4], zb[4], 'coral', zorder=-10)

#some shading
ax.fill_between(dist[0], zs[0], zb[0], color='royalblue', zorder=-10, alpha=0.1)
ax.fill_between(dist[0], bed[0], np.ones_like(bed[0])-3000, color='tan', zorder=-10, alpha=0.7)

#some plot params
ax.set_ylim(-2000,750)
ax.set_xlim(np.min(dist),np.max(dist))
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Elevation (m)')
ax.legend()


#%% Show the flowline on the map
fig2, ax1 = Plot_Antarctica2(GL=True, basemap = 'light', continental_shelf=0.0, extent=map_extent)
ax1.scatter(x[0],y[0], c=bed[0]-zb[0])

