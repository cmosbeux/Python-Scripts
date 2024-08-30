#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:33:47 2022

@author: cmosbeux

Description: Functions for vtu comparisons
"""

import pyvista
import numpy as np
import copy
from scipy.interpolate import griddata
from shapely.geometry import Polygon
import os

import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
import ParticleTracking
from scipy.interpolate import griddata
import time



#%% Make some masking
def out_mask(A):
    minx, maxx, miny,maxy = A[0].min(), A[0].max(), A[1].min(), A[1].max()
    ext_bound = Polygon([(minx,miny), (minx,maxy), (maxx,miny), (maxx,maxy)])
    A = A.transpose()
    in_bound = Polygon(A) 
    return ext_bound


def concat(*arrs) -> np.ndarray:
    return np.concatenate(tuple(map(np.asarray, arrs)))

def insert_at(outer_arr, arr, n) -> np.ndarray:
    outer_arr = np.asarray(outer_arr)
    prev, post = np.split(outer_arr, (n,))
    return concat(prev, arr, post)


def cross2d(x1, y1, x2, y2):
    return x1*y2-x2*y1


def is_clockwise(x1, y1, x2, y2):
    cp = cross2d(x1, y1, x2, y2)
    return cp < 0 if cp != 0 else None


def fill_outside(x, y, ll, ur, counter_clockwise=None):
    """
    Creates a polygon where x and y form a crevice of an outer
    rectangle with lower left and upper right corners `ll` and `ur`
    respectively. If `counter_clockwise` is `None` then the orientation
    of the outer polygon will be guessed to be the opposite of the
    inner connecting points.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    xmin, ymin = ll
    xmax, ymax = ur
    xmin, ymin = min(xmin, min(x)), min(ymin, min(y))
    xmax, ymax = max(xmax, max(x)), max(ymax, max(y))
    corners = np.array([
        [xmin, ymin],
        [xmin, ymax],
        [xmax, ymax],
        [xmax, ymin],
        [xmin, ymin],
    ])
    lower_left = corners[0]
    # Get closest point to splicing corner
    x_off, y_off = x-lower_left[0], y-lower_left[1]
    closest_n = (x_off**2+y_off**2).argmin()
    # Guess orientation
    p = [x_off[closest_n], y_off[closest_n]]
    try:
        pn = [x_off[closest_n+1], y_off[closest_n+1]]
    except IndexError:
        # wrap around if we're at the end of the array
        pn = [x_off[0], y_off[0]]
    if counter_clockwise is None:
        counter_clockwise = not is_clockwise(*p, *pn)
    corners = corners[::-1] if counter_clockwise else corners
    # Join the arrays
    corners = concat(np.array([[x[closest_n], y[closest_n]]]), corners)
    xs, ys = np.transpose(corners)
    return insert_at(x, xs, closest_n), insert_at(y, ys, closest_n)


def get_amundsen_contour():
    print('loading contour files...')
    Contour_dir = '../Contour_Amundsen/'
    A1=np.loadtxt(Contour_dir+'Contour_2.txt')
    A2=np.loadtxt(Contour_dir+'Contour_1.txt')
    A1 = A1.transpose()
    A2 = A2.transpose()
    A=np.concatenate((A1,A2), axis = 1)
    return A
    

def maskout(ax, A):
    minx, maxx, miny,maxy = A[0].min(), A[0].max(), A[1].min(), A[1].max()
    x_fill, y_fill = fill_outside(A[0], A[1], (minx,miny), (maxx,maxy))
    return ax.fill(x_fill, y_fill, color = 'w', zorder=10)


def outline_to_mask(line, x, y):
    """Create mask from outline contour

    Parameters
    ----------
    line: array-like (N, 2)
    x, y: 1-D grid coordinates (input for meshgrid)

    Returns
    -------
    mask : 2-D boolean array (True inside)
    """
    import matplotlib.path as mplp
    mpath = mplp.Path(line)
    X, Y = np.meshgrid(x, y)
    points = np.array((X.flatten(), Y.flatten())).T
    mask = mpath.contains_points(points).reshape(X.shape)
    return mask


#%%
def vtu2grid(filename, map_extent, dx=5e2, mask=True):
    print('Opening %s' % filename)
    reader = pyvista.get_reader(filename)
    mesh = reader.read()

    #get boundaries
    boundaryID = mesh.cell_data['GeometryIds']
    boundaryIndex = np.where(boundaryID == 104) #surface
    data = mesh.extract_cells(boundaryIndex)

    #surface data
    coords = data.points
    dZsdt = data.point_data['dsdt 3']
    v = data.point_data['velocity'].T
    vx, vy = v[0], v[1]
    velocity = np.sqrt(vx**2+vy**2)
    vobs = np.sqrt(data.point_data['vx']**2+data.point_data['vx']**2)

    #plot gmask (on base)
    boundaryIndex = np.where(boundaryID == 103)  #bed
    data = mesh.extract_cells(boundaryIndex)
    gm_stokes = data.point_data['groundedmask']

    melt = data.point_data['melt']
    #melt_ssa= data.point_data['meltssa']
    dZbdt = data.point_data['dsdt 3']
    H = data.point_data['thickness']
    vobs = np.sqrt(data.point_data['vx']**2+data.point_data['vy']**2)
    
    #make the dhdt computation
    dhdt_stokes = dZsdt - dZbdt
    
    #gl metrics
    x_gl, y_gl, v_gl, h_gl = coords[:,0][gm_stokes==0], coords[:,1][gm_stokes==0], velocity[gm_stokes==0], H[gm_stokes==0]

    #interpolate on a grid
    grid_x, grid_y = np.mgrid[map_extent[0]:map_extent[1]:dx, map_extent[2]:map_extent[3]:dx]
    
    if mask:
        if os.path.isfile('./Amundsen_mask.npy'):
            print('\t mask found...')
            mask1 = np.load('./Amundsen_mask.npy')
        else:
            print('\t making mask...')
            A = get_amundsen_contour()     
            mask1 = outline_to_mask(A.T, grid_x.T[0], grid_y[0])
            mask1=~mask1[::-1,:]
            np.save('./Amundsen_mask.npy', mask1)
            
    
    #gridding on regular grid
    print('\t Starting interpolation at resolution %0.2f \n' % (dx))
    melt_grid = griddata(coords[:,:2], melt, (grid_x, grid_y))[:,::-1].T
    #meltssa_grid = griddata(coords[:,:2], melt_ssa, (grid_x, grid_y))[:,::-1].T
    gm_grid = griddata(coords[:,:2], gm_stokes, (grid_x, grid_y))[:,::-1].T
    dhdt_grid = griddata(coords[:,:2], dhdt_stokes, (grid_x, grid_y))[:,::-1].T
    velocity_grid = griddata(coords[:,:2], velocity, (grid_x, grid_y))[:,::-1].T
    vobs_grid = griddata(coords[:,:2], vobs, (grid_x, grid_y))[:,::-1].T
    H_grid = griddata(coords[:,:2], H, (grid_x, grid_y))[:,::-1].T
    
    gm_grid = np.ma.masked_array(gm_grid, mask1)
    dhdt_grid = np.ma.masked_array(dhdt_grid, mask1)
    velocity_grid = np.ma.masked_array(velocity_grid, mask1)
    vobs_grid = np.ma.masked_array(vobs_grid,mask1)
    melt_grid = np.ma.masked_array(melt_grid,mask1)
    H_grid = np.ma.masked_array(H_grid, mask1)
    #some trick to make a grounded mask from SSA
    # gmssa_grid = copy.copy(meltssa_grid)
    # gmssa_grid[gmssa_grid>0] = -1 #ungrounded
    # gmssa_grid[gmssa_grid==0] = 1  #grounded
    # gmssa_grid[gmssa_grid<-1] = np.nan #out of bounds
    return gm_grid, dhdt_grid, melt_grid, velocity_grid, vobs_grid, H_grid #, gmssa_grid, meltssa_grid


def vtu2section(filename, map_extent, section, length=150,dx=5e2, mask=True , saveflowline = False ,fixed_flowline = False):
    '''This function reads a vtu, interpolates it on a regular grid and plot a flowline starting from a given point
    The user can save the flow line if desired.
    
    2 types of entry:
        - section = [x0, y0] --> starting point of a flowline that will be computed (or already exist if you saved it before)
        - section = [x0,y0,x1,y1] --> starting and ending point of a section.
    '''
    
    print('Opening %s' % filename)
    reader = pyvista.get_reader(filename)
    mesh = reader.read()

    #get boundaries
    boundaryID = mesh.cell_data['GeometryIds']
    boundaryIndex = np.where(boundaryID == 104) #surface
    data = mesh.extract_cells(boundaryIndex)

    #surface data
    coords = data.points
    dZsdt = data.point_data['dsdt 3']
    zs = data.point_data['fs upper']
    v = data.point_data['velocity'].T
    vx, vy = v[0], v[1]
    velocity = np.sqrt(vx**2+vy**2)
    vobs = np.sqrt(data.point_data['vx']**2+data.point_data['vx']**2)

    #plot gmask (on base)
    boundaryIndex = np.where(boundaryID == 103)  #bed
    data = mesh.extract_cells(boundaryIndex)
    gm_stokes = data.point_data['groundedmask']
    melt = data.point_data['melt']
    dZbdt = data.point_data['dsdt 3']
    H = data.point_data['thickness']
    
    zb = data.point_data['fs lower']
    bed = data.point_data['bed']
    
    #make the dhdt computation
    dhdt_stokes = dZsdt - dZbdt
    
    grid_x, grid_y = np.mgrid[map_extent[0]:map_extent[2]:dx, map_extent[1]:map_extent[3]:dx]
    
    #mask
    if mask:
        if os.path.isfile('./Amundsen_mask.npy'):
            print('\t mask found...')
            mask1 = np.load('./Amundsen_mask.npy')
            mask1 = mask1[::-1,:]
        else:
            print('\t making mask...')
            A = get_amundsen_contour()     
            mask1 = outline_to_mask(A.T, grid_x.T[0], grid_y[0])
            mask1=~mask1[::-1,:]
            np.save('./Amundsen_mask.npy', mask1)
            
    
    #gridding on regular grid
    print('\t Starting interpolation at resolution %0.2f \n' % (dx))
    zs_grid = griddata(coords[:,:2], zs, (grid_x, grid_y)).T
    zb_grid = griddata(coords[:,:2], zb, (grid_x, grid_y)).T
    bed_grid = griddata(coords[:,:2], bed, (grid_x, grid_y)).T
    velocity_grid = griddata(coords[:,:2], velocity, (grid_x, grid_y)).T
    vobs_grid = griddata(coords[:,:2], vobs, (grid_x, grid_y)).T
    vx_grid = griddata(coords[:,:2], vx, (grid_x, grid_y)).T
    vy_grid =  griddata(coords[:,:2], vy, (grid_x, grid_y)).T
    
    z = np.sqrt((grid_x/1e6)**2 + (grid_y/1e6)**2) + np.sin((grid_x/1e6)**2 + (grid_y/1e6)**2)

    # Coordinates of the line we'd like to sample along
    if len(section) == 4:
        x0,y0,x1,y1 = section 
        line = [(x0,y0), (x1,y1)]
        
        # Convert the line to pixel/index coordinates
        x_world, y_world = np.array(list(zip(*line))) 
        col = z.shape[1] * (x_world - grid_x.min()) / grid_x.ptp()    
        row = z.shape[0] * (y_world - grid_y.min()) / grid_y.ptp()
        
        # Interpolate the line at "num" points...
        num = 1000
        row, col = [np.linspace(item[0], item[1], num) for item in [row, col]]
        
        # Extract the values along the line, using cubic interpolation
        zs = scipy.ndimage.map_coordinates(np.nan_to_num(zs_grid), np.vstack((row, col)))
        zb = scipy.ndimage.map_coordinates(np.nan_to_num(zb_grid), np.vstack((row, col)))
        bed = scipy.ndimage.map_coordinates(np.nan_to_num(bed_grid), np.vstack((row, col)))
        x = scipy.ndimage.map_coordinates(np.nan_to_num(grid_x), np.vstack((row, col)))
        y = scipy.ndimage.map_coordinates(np.nan_to_num(grid_y), np.vstack((row, col)))
        velocity = scipy.ndimage.map_coordinates(np.nan_to_num(vobs_grid), np.vstack((row, col)))
        
    #Coordinates of starting point for downstream flowline
    elif len(section) == 2:       
        if fixed_flowline:
            x = np.loadtxt('flowline_x.txt')
            y = np.loadtxt('flowline_y.txt')
        else:
            p = section
            V = np.array([vx_grid, vy_grid])
            t_step = 2
            n_iter = round(length/2)
            x, y = ParticleTracking.VF_trajectory(p, grid_x.T, grid_y.T, V, t_step, n_iter)
            np.savetxt('flowline_x.txt', x)
            np.savetxt('flowline_y.txt', y)
        
        #flattening
        grid_x = grid_x.flatten()
        grid_y = grid_y.flatten()
        zs_grid= zs_grid.T.flatten()
        zb_grid = zb_grid.T.flatten()
        bed_grid = bed_grid.T.flatten()
        vobs_grid = vobs_grid.T.flatten()
        
        #reduce region to flowline bounds 
        cond_x = np.logical_and(grid_x<x.max(), grid_x>x.min()) 
        cond_y = np.logical_and(grid_y<y.max(), grid_y>y.min()) 
        cond = np.logical_and(cond_x,cond_y)
        
        grid_x = grid_x[cond]
        grid_y = grid_y[cond]
        zs_grid = zs_grid[cond]
        zb_grid = zb_grid[cond]
        bed_grid = bed_grid[cond]
        vobs_grid = vobs_grid[cond]
        
        #check plot
        # plt.scatter(grid_x.flatten(), grid_y.flatten(), c=zs_grid.T.flatten())
        # plt.plot(x,y,zorder = 50)
        
        #interpolate value 
        points = list(zip(grid_x, grid_y))
        # zs_flat = zs_grid.T.flatten()
        # zb_flat = zb_grid.T.flatten()
        # bed_flat = bed_grid.T.flatten()
        
        
        xflow, yflow, zs, zb, bed = [],[],[], [],[]
        vobs = []
        k=0
        print('\t interpolating data over the flowline')
        for xyi in zip(x,y):  
            k+=1
            
            #if we want mask and dist to vtu point too high
            cond1 = mask
            cond2 = np.min((xyi[0] - coords[:,0])**2 +  (xyi[1] - coords[:,1])**2)**0.5
            if cond2<10e3:
                xflow.append(xyi[0])
                yflow.append(xyi[1])
                zs.append(griddata(points, zs_grid, xyi, method='nearest'))
                zb.append(griddata(points, zb_grid, xyi, method='nearest'))
                bed.append(griddata(points, bed_grid, xyi,  method='nearest'))
                vobs.append(griddata(points, vobs_grid, xyi, method='nearest'))
            else:
                print('removing point from the flowline, %d m away...' % cond2)
        
    else:
        print('Wrong section/point input. It must be a [x0,y0,x1,y1] section or a [x0,y0] point...')
    
    # # Plot...
    # fig, axes = plt.subplots(nrows=2)
    # axes[0].pcolormesh(grid_x, grid_y, bed_grid.T,vmax=3000)
    # axes[0].plot(x_world, y_world, 'ro-')
    # axes[0].axis('image')
    
    # axes[1].plot(zs, c='cyan')
    # axes[1].plot(zb, c='blue')
    # axes[1].plot(bed, 'k')
    
    #save the flow line?
    if saveflowline: 
        np.savetxt('Flowline/x.txt', xflow)
        np.savetxt('Flowline/y.txt', yflow)
        
        #distance to top
        xflow = np.array(xflow)
        yflow = np.array(yflow)
        zs = np.array(zs)
        zb = np.array(zb)
        bed = np.array(bed)
        vobs = np.array(vobs)
        dist = ((xflow[0]-xflow)**2+(yflow[0]-yflow)**2)**0.5
        
        np.savetxt('Flowline/zs.txt', np.array([dist,zs]).T)
        np.savetxt('Flowline/zb.txt', np.array([dist,zb]).T)
        np.savetxt('Flowline/bed.txt', np.array([dist,bed]).T)
        np.savetxt('Flowline/vobs.txt',np.array([dist,vobs]).T)
    
    return xflow, yflow, zs, zb, bed


