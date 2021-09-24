# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 10:51:41 2015

@author: cmosbeux

Description:
------------

Functions to apply map projections as defined by the U.S. Geological Survey" 1984 pp 162-164.
The psxy2ll was coded by Craig Stewart, British Antarctic Survey 7/8/02 and adapted to python by C. Mosbeux (2015).
"""

import numpy as np

def psxy2ll(x,y):     
    '''
    Converts polarstereographic coordinates into WGS 84 Lat lon, 
    assuming a projection plane latitude of 71 South
    and a centre Meridian of 0 West.

    Reference: Snyder, J. P. 
    "Map projections Used by the U.S. Geological Survey" 1984 pp 162-164
    Coded by Craig Stewart, British Antarctic Survey 7/8/02
    adapted to python by C. Mosbeux (2015)
    Vectorised by Adrian Jenkins, British Antarctic Survey 22/7/02
    
    Define semi-major axis and flatteing of ellipsoid and
    standard parallel and central meridian for projection
    '''
    
    
    sin, cos, arctan, tan, arctan2 = np.sin, np.cos, np.arctan, np.tan, np.arctan2
    pi, sqrt = np.pi, np.sqrt
    
    a = 6378137                #WGS-84 major axis
    f = 1/298.257223563        #WGS-84 flattening
    
    Sp = -71                   #ADD projection plane latitude for southern polar stereographic (latitude of true scale)
    Cm = 0                     #Centre meridian
    phi_c = -Sp*pi/180         #Using Snyders nomenclature
    
    e = sqrt(2*f-f**2)
    
    tol = 1e-12; # set the latitude tolerance for iteration
    tc = tan(pi/4-phi_c/2)/((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(e/2) # eqn 13-9
    mc = cos(phi_c)/(1-e**2*(sin(phi_c))**2)**0.5 # eqn 12-15
    rho = sqrt(x**2 + y**2) # eqn 16-18
    t = rho*tc/(a*mc) # eqn 17-40
    lamda = Cm*pi/180 + arctan2(-x,y)       # eqn 16-16
    #lamda = -Cm*pi/180 + atan2(-x,-y);     # eqn 16-16
    
    # Have to iterate to find lat
    # First go use...
    phi1 = pi/2-2*arctan(t);
    trial_phi = pi/2 - 2*arctan(t*((1-e*sin(phi1))/(1+e*sin(phi1)))**(e/2))
    the_change = 2*tol # arbitrary number larger than 
    loop = True
    while loop == True:
        old_phi = trial_phi
        trial_phi = pi/2 - 2*arctan((t*((1-e*sin(old_phi))/(1+e*sin(old_phi)))**(e/2))) # eqn 7-9 (more or less)
        the_change = trial_phi-old_phi
        #print size(the_change), size(tol)
        if ((the_change) < tol):
            loop = False
            
    phi = trial_phi;
    # turn lat and lon into degrees
    lat = -phi*180/pi
    lon = -lamda*180/pi
    
    return lat, lon


def ll2psxy(lat,lon,Sp=-90.0, Cm=0, k0=0.97276901289):

    sin, cos, tan = np.sin, np.cos, np.tan
    pi, sqrt = np.pi, np.sqrt    
        
    """give lat, lon, Sp (latitude of true scale) and Cm (Centre meridian)"""
    phi_c = -Sp*pi/180       #Using Snyders nomenclature
    lamda0 = Cm*pi/180
    a=6378137.0              #WGS-84 major axis
    f=1.0/298.257223563     #WGS-84 flattening
    e = sqrt(2*f-f**2)
    
    phi=-lat*pi/180
    lamda=-lon*pi/180
    
    t  = tan(pi/4-phi/2)/((1-e*sin(phi))/(1+e*sin(phi)))**(e/2)
    tc = tan(pi/4-phi_c/2)/((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(e/2)
    mc = cos(phi_c)/(1-e**2*(sin(phi_c))**2)**0.5 
    
    #compute rho
    if Sp==-90.0:
        rho = 2*a*k0*t/((1+e)**(1+e)*(1-e)**(1-e))**0.5 
    else:
        rho= a*mc*t/tc*k0
   
    x=-rho*sin(lamda-lamda0)
    y=+rho*cos(-lamda-lamda0)
    
    return x,y

