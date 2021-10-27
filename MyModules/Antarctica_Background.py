#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 08:42:22 2019

@author: cmosbeux

Description:
------------
This file contains various functions to plot Antarctic maps
with grounding lines, observed velocities, continental shelf line, a scale,...

The main function is Plot_Antarctica that can be called as follow:
ax = Plot_Antarctica(args)

You can then apply some extra work to the ax by calling it ax[i], where i is
the axe index (e.g., i=0 if only one subplot, i=0,1 if two subplots, ...).
"""

import numpy as np
import matplotlib.pyplot as plt
import shapefile
from imageio import imread
from mpl_toolkits.axes_grid1 import AxesGrid
import csv
import os

#get current directory
dirname = os.path.dirname(__file__)
#%%

def plot_GL(ax, color='black', lw=0.5, icefront=True, precision=1,zorder=1000000):
    path = 'MyModuleData/GL/scripps_antarctica_polygons_v1.shp'  
    src_file = os.path.join(dirname, path)
    shp = shapefile.Reader(src_file)
    k=0
    for shape in shp.shapeRecords():
        if 'shelf' not in shape.record[1]:
            xy = [i for i in shape.shape.points[:]]
            x,y = zip(*[(j[0],j[1]) for j in xy])
            ax.plot(x,y, color=color, linewidth=lw, zorder=zorder)
        if 'shelf' in shape.record[1] and icefront:
            xy = [i for i in shape.shape.points[:]]
            xs,ys = zip(*[(j[0],j[1]) for j in xy])
            xxs, yys = np.zeros(len(xs))*np.nan, np.zeros(len(ys))*np.nan
            for i in np.arange(0,len(xs)-1):
                d2 = (xs[i+1]-xs[i])**2+(ys[i+1]-ys[i])**2
                if d2<1e8:
                    xxs[i] = xs[i]
                    yys[i] = ys[i]
            ax.plot(xxs,yys, color=color, linewidth=lw, zorder=zorder)
                    

def save_GL2021(color='black', lw=0.5, icefront=True, precision=1,zorder=1000000):
    shape_name=os.path.expanduser('MyModuleData/GL/scripps_antarctica_polygons_v1.shp')    
    shp = shapefile.Reader(shape_name)
    k=0
    x_save,y_save = [],[]
    save = []
    for shape in shp.shapeRecords():
        if 'shelf' not in shape.record[1]:
            xy = [i for i in shape.shape.points[:]]
            x,y = zip(*[(j[0],j[1]) for j in xy])
            x_save.append(x)
            y_save.append(y)
            save.append(zip(x,y))
        if 'shelf' in shape.record[1] and icefront:
            # xs=np.asarray(shape.shape.points[:]).transpose()[0][::precision]
            # ys=np.asarray(shape.shape.points[:]).transpose()[1][::precision]
            xy = [i for i in shape.shape.points[:]]
            xs,ys = zip(*[(j[0],j[1]) for j in xy])
            xxs, yys = np.zeros(len(xs))*np.nan, np.zeros(len(ys))*np.nan
            for i in np.arange(0,len(xs)-1):
                d2 = (xs[i+1]-xs[i])**2+(ys[i+1]-ys[i])**2
                if d2<1e8:
                    xxs[i] = xs[i]
                    yys[i] = ys[i]            
    return x_save,y_save
            
              
def plot_continental_shelf(ax, color='dimgrey', lw=0.5, precision=1,zorder=1000000):
    path = 'MyModuleData/BATHY/Bathymetry_contours/contour_1500m_simple.npy'
    src_file = os.path.join(dirname, path)
    A = np.load(src_file)
    x, y = A[0], A[1]
    ax.plot(x,y, color=color, lw=lw)


def plot_continental_shelf2(ax, color='darkgrey', lw=0.0, precision=1,zorder=1000000):
    path = 'MyModuleData/BATHY/Bathymetry_contours/contour_1500m_simple.npy'
    src_file = os.path.join(dirname, path)
    with open(src_file) as csvfile:
        reader = csv.reader(csvfile, delimiter = '\n')
        x, y = [], []
        for row in reader:
            xy = row[0].split(',')
            x.append(float(xy[0]))
            y.append(float(xy[1]))
    ax.scatter(x[::precision],y[::precision],c=color, s=0.5,linewidth=lw, zorder=zorder)
 
    
 
    
 
#%% 
        
def plot_front(ax, color='black', lw=1, zorder=1000000):
    try:
        path = 'MyModuleData/ICE_FRONT/ice_front_antarctica_xy.npy'
        src_file = os.path.join(dirname, path)
        x,y = np.load(src_file)
    except IOError:
        print('Error', 'numpy file does not exist... creating it!')
        try:
            from SouthPolar_Coordinates import ll2psxy
            xx, yy = [], []
            path = 'MyModuleData/ICE_FRONT/ice_front_antarctica_lat_lon.txt'
            src_file = os.path.join(dirname, path)
            with open(src_file) as f:
                for line in f:
                    if '#' in line or line == '\n':
                        continue
                    else:
                        d = line.split()
                        try:
                            lat, lon = float(d[0]), float(d[1]) 
                            x, y = ll2psxy(lat,lon)
                            xx.append(x), yy.append(y)
                        except ValueError:
                            None
                        except IndexError:
                            None
                xx, yy = np.asarray(xx), np.asarray(yy)
                A = np.array([xx,yy])
                path = 'MyModuleData/ICE_FRONT/ice_front_antarctica_xy'
                src_file = os.path.join(dirname, path)
                np.save(src_file, A) 
                
        except IOError:
            print('Error', 'source txt file does not exist... check sources and path.')

    ax.scatter(x,y, color = 'k', s=0.5, zorder=zorder)   
    return x,y

def velocity_background(ax, color='white'):         
    np.random.seed(0)
    path = 'MyModuleData/IMAGE/RIS_veloctiy_background.png'
    src_file = os.path.join(dirname, path)
    img = imread(src_file)
    left=-1246954.0
    right=1503412.63317
    top=273093.6875
    bottom=-1363363.5765199999
    ax.imshow(img ,extent=[left, right, bottom, top], zorder=0)
    
    
def basemap_velocity(ax, cmap = 'Greys_r', norm = None, mask_null = False,alpha=0.5,zorder=0):         
    np.random.seed(0)
    path = './MyModuleData/IMAGE/Antarctica_velocity_basemap_grey.png'
    src_file = os.path.join(dirname, path)
    img = imread(src_file)
    if mask_null:
        img = np.ma.masked_array(img, mask=img==255)
    left =   -2811000.0
    right =  +2899500.0
    top =    +2767500.0
    bottom = -2345500.0
    if norm == None:
        ax.imshow(img,  extent=[left, right, bottom, top], interpolation = 'bicubic', cmap = cmap,  alpha=alpha, zorder=zorder)
    else:
        ax.imshow(img,  extent=[left, right, bottom, top], interpolation = 'bicubic', cmap = cmap,  norm=norm, alpha=alpha, zorder=zorder)
        
        
def basemap_LIMA(ax, color='white', cmap = 'Greys_r'):  
    np.random.seed(0)
    path = './MyModuleData/IMAGE/Antarctica_MOSAIC_basemap_grey.png'
    src_file = os.path.join(dirname, path)       
    img = imread(src_file)
    left =   -2811000.0
    right =  +2899500.0
    top =    +2767500.0
    bottom = -2345500.0
    ax.imshow(img,  extent=[left, right, bottom, top], interpolation = 'bicubic', cmap = cmap,  alpha=1, zorder=0)
    

def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def cbar(fig, mappable, ticks=None, tickslabel=None, orientation='horizontal'):
    cbaxes = fig.add_axes([0.25, 0.15, 0.2, 0.02]) 
    cbar = fig.colorbar(mappable, cax=cbaxes, ticks=ticks, orientation=orientation)
    if tickslabel is not None:
        cbar.ax.set_xticklabels(tickslabel, fontsize=16)
    return  cbar


def Plot_Antarctica(nrows=1, ncols=1, GL=True, icefront=True, continental_shelf=0, precision=1, basemap=None, extent=[-3333000, 3333000, -3333000, 3333000], cbar=None, axes_pad=0.1, figsize=(20,20)):
    """Extent should be given: (xmin, xmax, ymin, ymax)
       a basemap ca be choosen: dark, light, or blue.
    """
    
    xmin, xmax, ymin, ymax = extent
    fig = plt.figure(figsize=figsize)
    #fig, ax= plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,8*nrows))
        
    if cbar == None:
        grid = AxesGrid(fig, 111,
            nrows_ncols=(nrows, ncols),
            axes_pad=axes_pad,
            cbar_mode=None,
            cbar_location='right',
            )  
    elif cbar == 'single':
        grid = AxesGrid(fig, 111,
            nrows_ncols=(nrows, ncols),
            axes_pad=axes_pad,
            cbar_mode='single',
            cbar_location='right',
            cbar_pad=0.05,
            cbar_size = "5%",
            )     
    elif cbar == 'each':
        grid = AxesGrid(fig, 111,  # similar to subplot(143)
            nrows_ncols=(nrows, ncols),
            axes_pad=0.75,
            label_mode="1",
            share_all=False,
            cbar_location="bottom",
            cbar_mode="each",
            cbar_size="4%",
            cbar_pad="2%",
            )
   
    for ax in grid:
        ax.set_xticks([]), ax.set_xticklabels([])
        ax.set_yticks([]), ax.set_yticklabels([])
        if GL:
            if isinstance(GL,int):
                plot_GL(ax, icefront=icefront, lw = GL, precision = precision)
            else:
                #default GL linewidth = 0.5
                plot_GL(ax, icefront=icefront, precision = precision)
            #plot continental shelf
            plot_continental_shelf(ax, lw=continental_shelf)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        if basemap == 'dark':
            basemap_velocity(ax, cmap='Greys', alpha = 1)
        elif basemap == 'light':
            basemap_velocity(ax, cmap='Greys_r')
        elif basemap == 'blue':
            basemap_velocity(ax, cmap='Blues_r', alpha = 1)
        elif basemap is not None and 'blue' in basemap:
            try:
                a = float(basemap[-3:])
            except ValueError:
                a = 1
                print('Incorrect alpha value for blue velocity map, use default alpha = 1 ')
            basemap_velocity(ax, cmap='Blues_r' ,alpha = a)
        else:
            continue
    return fig, grid


def Plot_Antarctica2(nrows=1, ncols=1, GL=True, icefront=True, continental_shelf=0, precision=1, basemap=None, extent=[-3333000, 3333000, -3333000, 3333000], cbar=None, axes_pad=0.1, figsize=(20,20)):
    """Extent should be given: (xmin, xmax, ymin, ymax)"""
    
    xmin, xmax, ymin, ymax = extent
    #fig = plt.figure(figsize=figsize)
    fig, ax= plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols,8*nrows))
        
    def plot_basemap(basemap, ax):
        if basemap == 'dark':
            out = basemap_velocity(ax, cmap='Greys', alpha = 1)
        elif basemap == 'light':
            out = basemap_velocity(ax, cmap='Greys_r')
        elif basemap == 'blue':
            out = basemap_velocity(ax, cmap='Blues_r' ,alpha = 1)
        elif basemap == 'lightblue':
            out = basemap_velocity(ax, cmap='Blues_r' ,alpha = 0.8)
        elif basemap==None:
            out = 0
        else:
            out = print( "Requested bedmap >%s< does not exist")
        return out
    
    #one single figure
    if nrows == 1 and ncols == 1:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        for tickx, ticky in zip(ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()):
            tickx.label.set_fontsize(8) 
            ticky.label.set_fontsize(8)             
        plot_GL(ax)    
        plot_basemap(basemap, ax)
     
    #raw or col figure
    elif (ncols+nrows) > 2 and min(ncols,nrows)==1:
        for j in np.arange(0,max(ncols,nrows)):
            ax[j].set_xticks([])
            ax[j].set_yticks([])
            ax[j].set_xlim(xmin,xmax)
            ax[j].set_ylim(ymin,ymax)
            for tickx, ticky in zip(ax[j].xaxis.get_major_ticks(),ax[j].yaxis.get_major_ticks()):
                        tickx.label.set_fontsize(8) 
                        ticky.label.set_fontsize(8)           
            plot_GL(ax[j])
            plot_basemap(basemap,ax[j])
                    
    #multiple column and line figure            
    elif nrows>1 and ncols>1:
        for j in np.arange(0,nrows):
            for i in np.arange(0,ncols):
                ax[i][j].set_xticks([])
                ax[i][j].set_yticks([])
                ax[i][j].set_xlim(xmin,xmax)
                ax[i][j].set_ylim(ymin,ymax)
                for tickx, ticky in zip(ax[i][j].xaxis.get_major_ticks(),ax[i][j].yaxis.get_major_ticks()):
                            tickx.label.set_fontsize(8) 
                            ticky.label.set_fontsize(8)           
                plot_GL(ax[i][j])
                plot_basemap(basemap, ax[i][j])
                            
    #fig.tight_layout()
    return fig, ax


def scale(ax, length=100):
    """ Plot a scale on a map (default scale is 100 km)"""
    #scale
    l = length*1e3
    ll = '%d km' % length
    x = ax.get_xbound()[0]
    y = ax.get_ybound()[0]
    y = y + abs(y)*0.05
    delta = abs(x)*0.05
    x = x+delta
    xx = [x, x+l]
    
    ax.plot(xx, [y,y], c = 'black', zorder = 1e12)
    ax.text(x+l/2, y+l/4, s=ll, ha='center', color='black', fontweight='medium', fontsize=11, zorder = 1e12)


class basin:
    """extent of different basins"""
    @staticmethod
    def RIS():
        #xmin, xmax, ymin, ymax
        return [-6.0e5,4.5e5,-1.4e6,-4e5]
    def PanAntarctic():
        #xmin, xmax, ymin, ymax
        return [-3.0e6,3e6,-3e6,3e6]
    def Amundsen():
        #xmin, xmax, ymin, ymax
        return [-2e6,-1e6,-9e5,1e5]
        
    
   
