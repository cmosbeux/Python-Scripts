#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 22:11:37 2017

@author: cmosbeux
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

LinL = np.loadtxt('0-1/Linear_L_0-1.csv', delimiter=',')

b3=LinL[:,2] # value of blue at sample n
b2=LinL[:,2] # value of blue at sample n
b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1

# setting up columns for list
g3=LinL[:,1]
g2=LinL[:,1]
g1=np.linspace(0,1,len(g2))

r3=LinL[:,0]
r2=LinL[:,0]
r1=np.linspace(0,1,len(r2))

# creating list
R=zip(r1,r2,r3)
G=zip(g1,g2,g3)
B=zip(b1,b2,b3)

# transposing list
RGB=zip(R,G,B)
rgb=zip(*RGB)
# print rgb

# creating dictionary
k=['red', 'green', 'blue']
LinearL=dict(zip(k,rgb)) # makes a dictionary from 2 lists


my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',LinearL)
#plt.pcolor(np.random.rand(10,10),cmap=my_cmap)
#plt.colorbar()