#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 20:08:00 2017

@author: cmosbeux
"""
#!/usr/bin/env python
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap


cdict = {'red':   ((0.0, 1.0, 1.0), 
                   (0.1, 1.0, 1.0),  # red 
                   (0.4, 1.0, 1.0),  # violet
                   (1.0, 0.0, 0.0)), # blue
         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),  # red
                   (0.4, 1.0, 1.0),  # violet
                   (1.0, 1.0, 0.0))  # blue
          }

def alpha_cmap(cmap, type='mid', width=10):
    """apply some transparency to a colormap. 
    type = mid : alpha from the center
         = min : alpha from min
         = max : alpha from max """ 
    n = cmap.N
    width  = width//2*2
    width = max(width,2)
    #get the shape of the alpha map
    if type == 'mid':
        side = np.ones((n-width)//2)
        center = np.linspace(0,1,width//2)
        alpha_array = np.concatenate([side, center[::-1], center, side])   
    elif type == 'min':
        minl = np.linspace(0,1,width) 
        maxl = np.ones((n-width))
        alpha_array = np.concatenate([minl, maxl])   
    elif type == 'max':
        minl = np.linspace(0,1,width) 
        maxl = np.ones((n-width))
        alpha_array = np.concatenate([maxl, minl])   
    else: 
        print('warning: incorrect alpha map style, choose min, max or mid.')
        
    my_cmap = cmap(np.arange(cmap.N)) # Get the colormap colors
    my_cmap[:,-1] = alpha_array # Set alpha
    my_cmap = ListedColormap(my_cmap) # Create new colormap
    
    return my_cmap


def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)
       



#%% a few cmaps I created
c = mcolors.ColorConverter().to_rgb


def rvb():
    rvb=make_colormap(    
    [c('crimson'), c('midnightblue'), 0.1, 
     c('midnightblue'), c('royalblue'),0.20, 
     c('royalblue'), c('deepskyblue'), 0.30,
     c('deepskyblue'), c('cyan'), 0.40,
     c('cyan'), c('white'), 0.50,  
     c('white'), c('limegreen'), 0.65, 
     c('limegreen'), c('khaki'), 0.75,
     c('khaki'), c('sienna'), 0.95,
     c('sienna'), c('saddlebrown'), 1.0,
     c('saddlebrown')]
    )
    return rvb

def rvb2():
    rvb=make_colormap(    
    [c('crimson'), c('midnightblue'), 0.05, 
     #c('blueviolet'), c('darkslateblue'), 0.15, 
     c('midnightblue'), c('royalblue'),0.15, 
     c('royalblue'), c('deepskyblue'), 0.25,
     c('deepskyblue'), c('cyan'), 0.35,
     c('cyan'), c('limegreen'), 0.40,  
     #c('deepskyblue'), c('white'), 0.6,  
     c('limegreen'), c('white'), 0.50, 
     c('white'), c('khaki'), 0.75,
     c('khaki'), c('sienna'), 0.95,
     c('sienna'), c('saddlebrown'), 1.0,
     c('saddlebrown')]
    )
    return rvb

def nic_edge():
    cmap=make_colormap(    
     [c('midnightblue'), c('blue'),0.10, 
     c('blue'), c('cyan'), 0.35,
     c('cyan'), c('white'), 0.50,   
     c('white'), c('yellow'), 0.60,
     c('yellow'), c('red'), 0.90,
     c('red'), c('darkred'), 1.0,
     c('darkred')]
    )
    return cmap  
   
def Blues():
    cmap=make_colormap(    
     [c('midnightblue'), c('blue'),0.20, 
     c('blue'), c('cyan'), 0.70,
     c('cyan'), c('white'), 1.0,  c('white') ] 
     )
    return cmap


def dark_rainbow():
    yel=(1.0, 0.90, 0.0)
    dred=(0.9,0.1,0.2)
    dblue=(0.1,0.098,0.6)
    cmap=make_colormap(    
     [dblue, c('royalblue'), 0.2,
     c('royalblue'), c('limegreen'), 0.5,
     c('limegreen'), yel, 0.75,
     yel, c('orange'), 0.80,
     c('orange'), dred, 1.0,   
     dred]
    )
    return cmap  

# def terrain_div_2(self):
#     c=[(151,103,0)( 151,103,0)( 151,104,0)( 151,105,0)( 151,105,0)( 151,106,0)( 151,107,0)( 151,107,0)( 152,108,0)( 152,109,0)( 152,109,0)( 152,110,0)( 152,111,0)( 152,111,0)( 152,112,0)( 152,113,0)( 152,113,0)( 153,114,0)( 153,115,0)( 153,115,0)( 153,116,0)( 153,117,0)( 153,117,0)( 153,118,0)( 153,118,0)( 153,119,0)( 154,120,0)( 154,120,0)( 154,121,0)( 154,122,0)( 154,122,0)( 154,123,0)( 154,124,0)( 154,124,0)( 154,125,0)( 155,126,0)( 155,126,0)( 155,127,0)( 155,128,0)( 155,128,0)( 155,129,0)( 155,130,0)( 155,130,0)( 155,131,0)( 155,132,0)( 156,132,0)( 156,133,0)( 156,134,0)( 156,134,0)( 156,135,0)( 156,136,0)( 156,136,0)( 156,137,0)( 156,137,0)( 157,138,0)( 157,139,0)( 157,139,0)( 157,140,0)( 157,141,0)( 157,141,0)( 157,142,0)( 157,143,0)( 157,143,0)( 158,144,0)( 157,145,0)( 156,145,0)( 156,146,0)( 155,147,0)( 154,148,0)( 153,148,0)( 152,149,0)( 151,150,0)( 150,151,0)( 149,151,0)( 148,152,0)( 147,153,0)( 146,154,0)( 145,154,0)( 145,155,0)( 144,156,0)( 143,157,0)( 142,157,0)( 141,158,0)( 140,159,0)( 139,160,0)( 138,160,0)( 137,161,0)( 136,162,0)( 135,163,0)( 134,163,0)( 133,164,0)( 133,165,0)( 132,166,0)( 131,166,0)( 130,167,0)( 129,168,0)( 128,169,0)( 127,169,0)( 126,170,0)( 125,171,0)( 124,172,0)( 123,172,0)( 122,173,0)( 121,174,0)( 121,175,0)( 120,176,0)( 119,176,0)( 118,177,0)( 117,178,0)( 116,179,0)( 115,179,0)( 114,180,0)( 113,181,0)( 112,182,0)( 111,182,0)( 110,183,0)( 110,184,0)( 109,185,0)( 108,185,0)( 107,186,0)( 106,187,0)( 105,188,0)( 104,188,0)( 103,189,0)( 102,190,0)( 101,191,0)( 100,191,0)( 99,192,0)( 100,193,2)( 103,194,6)( 105,195,10)( 108,196,14)( 110,197,18)( 112,198,22)( 115,199,26)( 117,200,30)( 120,201,34)( 122,202,38)( 125,203,42)( 127,204,46)( 130,205,50)( 132,206,54)( 135,207,58)( 137,208,62)( 140,209,66)( 142,210,70)( 144,211,74)( 147,212,78)( 149,213,82)( 152,214,86)( 154,215,90)( 157,216,94)( 159,217,98)( 162,218,102)( 164,219,106)( 167,220,110)( 169,221,114)( 171,222,118)( 174,223,122)( 176,224,126)( 179,225,130)( 181,226,134)( 184,227,138)( 186,228,142)( 189,229,146)( 191,230,150)( 194,231,154)( 196,232,158)( 199,233,162)( 201,234,167)( 203,235,171)( 206,236,175)( 208,237,179)( 211,238,183)( 213,239,187)( 216,240,191)( 218,241,195)( 221,242,199)( 223,243,203)( 226,244,207)( 228,245,211)( 231,246,215)( 233,247,219)( 235,248,223)( 238,249,227)( 240,250,231)( 243,251,235)( 245,252,239)( 248,253,243)( 250,254,247)( 253,255,251)( 255,255,255)( 253,253,255)( 249,249,255)( 245,245,254)( 241,241,254)( 237,237,253)( 233,233,253)( 229,229,252)( 225,225,252)( 221,221,251)( 217,217,250)( 213,213,250)( 209,209,249)( 205,205,249)( 201,201,248)( 197,197,248)( 193,193,247)( 189,189,247)( 185,185,246)( 181,181,246)( 177,177,245)( 173,173,244)( 169,169,244)( 164,164,243)( 160,160,243)( 156,156,242)( 152,152,242)( 148,148,241)( 144,144,241)( 140,140,240)( 136,136,240)( 132,132,239)( 128,128,238)( 124,124,238)( 120,120,237)( 116,116,237)( 112,112,236)( 108,108,236)( 104,104,235)( 100,100,235)( 96,96,234)( 92,92,234)( 88,88,233)( 84,84,232)( 80,80,232)( 76,76,231)( 72,72,231)( 68,68,230)( 64,64,230)( 60,60,229)( 56,56,229)( 52,52,228)( 48,48,228)( 44,44,227)( 40,40,226)( 36,36,226)( 32,32,225)( 28,28,225)( 24,24,224)( 20,20,224)( 16,16,223)( 12,12,223)( 8,8,222)( 4,4,222)( 0,0,221)]
#     c = np.array(c)
#     cm = mcolors.ListedColormap(c/255.)
#     return cm