#Lecture files
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv, math
import re
import string
import os
from os import listdir
from os.path import isfile, isdir, join
from matplotlib.collections import LineCollection

# files reading and creation of datas
#-------------------------------------

"""lecture of contour file : grounding line or catchement contour"""
def lecture_gl(arg):
    f = open(arg, 'rb') 
    reader = csv.reader(f, dialect = 'excel', delimiter=' ', skipinitialspace=True) 
    long_gl = []
    lat_gl = []
    for row in reader:
        long_gl.append(row[1])
        lat_gl.append(row[0])
        
    long_gl=np.asarray(map(float, long_gl[0:np.size(long_gl)]))
    lat_gl=np.asarray(map(float, lat_gl[0:np.size(lat_gl)]))
        
    return long_gl, lat_gl

def lecture_contour(delimiter,*args):
    x_contour = []
    y_contour = []
    for contour_file in args:
        f = open(contour_file, 'rt') 
        reader = csv.reader(f, dialect = 'excel', delimiter=delimiter) 
        for row in reader:
            x_contour.append(float(row[0]))
            y_contour.append(float(row[1]))
    #x_contour=np.asarray(map(float, x_contour[0:np.size(x_contour)]))
    #y_contour=np.asarray(map(float, y_contour[0:np.size(y_contour)]))  
    x_contour, y_contour = np.asarray(x_contour), np.asarray(y_contour)
    return x_contour, y_contour
    
    
# """----------------------------------------------------------------------------
# This function allows read the Vaughan Data for Antarctica and create .npy files
# which are lighter to load for plotting. A .csv file is also created.

# files to load are in 'CONVERT folder.

# """
# def lecture_data():
        
#     mypath = '/Users/cmosbeux/Documents/Python/East_Antarcitca_Project/Antarctica_Bedrock_Data'        
#     f_name = 'Mission-data-available-for-release.txt' 
    
#     f = open(mypath+'/'+f_name, 'rb')
    
#     header = f.readline().split(',') #read first line and split it
#     reader = csv.reader(f, delimiter=',') #read other lines from file
    
#     i=0
#     print('number of data field:')
#     while i<np.size(header)-1: #loop on header number
#         data = [] 
#         print(i)
#         for line in reader:
#             data.append(line[i])
#         exec("%s = %s" % (header[i], data)) #put header name of data     
#         i=i+1  
#         f.seek(0) #restart the lecture at the file beginning
#         a=f.readline() #pass header
        
#     i=0
#     n_max = np.size(Ice_thickness)
#     while i<n_max:
        
#         if Ice_thickness[i]=='':
#             Ice_thickness[i]='nan'
#         elif Ice_thickness[i]=='-9999':
#             Ice_thickness[i]='nan' 
               
#         if Surface_elev[i]=='':
#             Surface_elev[i]='nan'
#         elif Surface_elev[i]=='-9999':
#             Surface_elev[i]='nan'

#         if Bed_elev[i]=='':
#             Bed_elev[i]='nan'
#         elif Bed_elev[i]=='-9999':
#             Bed_elev[i]='nan'                                                    
       
#         i=i+1
        
#     #change list in array    
#     print('changing list to numpy array...')
#     longitude=np.asarray(map(float, Longitude[0:np.size(Longitude)]))
#     print('    longitude ok')
#     latitude=np.asarray(map(float, Latitude[0:np.size(Latitude)]))
#     print('    latitude ok')
#     thickness = np.asarray(map(float,Ice_thickness[0:np.size(Ice_thickness)]))
#     print('    thickness ok')
#     bedrock = np.asarray(map(float,Bed_elev[0:np.size(Bed_elev)]))
#     print('    bedrock ok')
#     surface = np.asarray(map(float,Surface_elev[0:np.size(Surface_elev)]))
#     print('    surface ok')
    
#     MATRIX = np.array([longitude, latitude, thickness, bedrock, surface])
#     MATRIX = np.transpose(MATRIX)
#     #np.savetxt('matrix.txt', MATRIX, delimiter=' ', newline='\n' ) 
#     np.save('data_Vaughan', MATRIX)
               
#     return longitude, latitude, thickness, bedrock, surface
    
    
    #Plot the density map using nearest-neighbor interpolation
    #plt.pcolormesh(X,Y,Z)
    #plt.show() 
    
"""----------------------------------------------------------------------------
This function allows read the Ice_Bridge_UTIG data for Antarctica and create 
.npy files which are lighter to load for plotting. A .csv file is also created.

input args* = Ice_Bridge_Operation/20XX_AN_UTIG
"""
   
def lecture_airborn_ice_bridge(*args):
    longitude = []
    latitude = []
    thickness = []
    bedrock = []
    surface = []
    #folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 

    with open('Ice_Bridge_UTIG.csv', 'w') as csvfile: 
        for mypath in args:
            print('Main Directory:', mypath)
            folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 
            for x in folders:
                print('    start reading directory: ', x)
                onlyfiles = [ f for f in listdir(mypath+'/'+x) if isfile(join(mypath+'/'+x,f)) ] 
                for f in onlyfiles :  
                    #if re.search('^ASB', f):
                    #print '        File :', f                                        
                    with open(mypath+'/'+x+'/'+f, 'r') as f: 
                        reader = csv.reader(f, dialect = 'excel', delimiter=' ')
                        for line in reader:
                            #for csv
                            csv_writer = csv.writer(csvfile, delimiter=',')  
                            csv_writer = csv_writer.writerow(line)  
                            #for numpy
                            longitude.append(line[3]) 
                            latitude.append(line[4])
                            thickness.append(line[5])
                            bedrock.append(line[7])
                            surface.append(line[8])
                                                                                                                                                    
    #change list in array    
    print('changing list to numpy array...')
    longitude=np.asarray(map(float, longitude[0:np.size(longitude)]))
    print('    longitude ok')
    latitude=np.asarray(map(float, latitude[0:np.size(latitude)]))
    print('    latitude ok')
    thickness = np.asarray(map(float,thickness[0:np.size(thickness)]))
    print('    thickness ok')
    bedrock = np.asarray(map(float,bedrock[0:np.size(bedrock)]))
    print('    bedrock ok')
    surface = np.asarray(map(float,surface[0:np.size(surface)]))
    print('    surface ok')
    
    MATRIX = np.array([longitude, latitude, thickness, bedrock, surface])
    MATRIX = np.transpose(MATRIX)
    #np.savetxt('matrix.txt', MATRIX, delimiter=' ', newline='\n' ) 
    np.save('data_IceBridge_UTIG', MATRIX)
               
    return longitude, latitude, thickness, bedrock, surface
    
def lecture_airborn_duncan(*args):
    longitude = []
    latitude = []
    thickness = []
    bedrock = []
    surface = []
    year=[]
    #folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 

    with open('Ice_Bridge_UTIG_Supplementary_date.csv', 'w') as csvfile: 
        for mypath in args:
            print('Main Directory:', mypath)
            folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 
            for x in folders:
                print('    start reading directory: ', x)
                onlyfiles = [ f for f in listdir(mypath+'/'+x) if isfile(join(mypath+'/'+x,f)) ] 
                for f in onlyfiles :
                    #print f
                    if re.search('.txt$', f):
                        print('          --> ', f)
                        with open(mypath+'/'+x+'/'+f, 'r') as f: 
                            reader = csv.reader(f, dialect = 'excel', delimiter=' ')
                            for line in reader:
                                #print line
                                if line[0]=='#': #re.search('^#', line):
                                    pass
                                else:                                
                                    #for csv
                                    csv_writer = csv.writer(csvfile, delimiter=',')  
                                    csv_writer = csv_writer.writerow(line)  
                                    #for numpy
                                    year.append(line[0])
                                    longitude.append(line[3]) 
                                    latitude.append(line[4])
                                    thickness.append(line[5])
                                    bedrock.append(line[7])
                                    surface.append(line[8])
                                                                                                                                                    
    #change list in array                                                                                                                                          print 'changing list to numpy array...'
    year=np.asarray(map(float, year[0:np.size(year)]))
    print('changing list to numpy array...')
    longitude=np.asarray(map(float, longitude[0:np.size(longitude)]))
    print('    longitude ok')
    latitude=np.asarray(map(float, latitude[0:np.size(latitude)]))
    print('    latitude ok')
    thickness = np.asarray(map(float,thickness[0:np.size(thickness)]))
    print('    thickness ok')
    bedrock = np.asarray(map(float,bedrock[0:np.size(bedrock)]))
    print('    bedrock ok')
    surface = np.asarray(map(float,surface[0:np.size(surface)]))
    print('    surface ok')
    
    MATRIX = np.array([longitude, latitude, thickness, bedrock, surface, year])
    MATRIX = np.transpose(MATRIX)
    #np.savetxt('matrix.txt', MATRIX, delimiter=' ', newline='\n' ) 
    np.save('data_IceBridge_UTIG_Supplementary_date', MATRIX)
               
    return year, longitude, latitude, thickness, bedrock, surface
               
def lecture_airborn_ice_bridge_bis(*args):#(mypath1, mypath2):
    longitude = []
    latitude = []
    thickness = []
    bedrock = []
    surface = []
    #folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 
    filename = 'IceBridge_IRMCR_3'
    with open(filename+'.csv', 'w') as csvfile: 
        for mypath in args:
            print('Main Directory:', mypath)
            folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 
            for x in folders:
                if re.search('^2013', x) or re.search('^2014', x):
                    print('    start reading directory: ', x)
                    onlyfiles = [ f for f in listdir(mypath+'/'+x) if isfile(join(mypath+'/'+x,f)) ] 
                    for f in onlyfiles :
                        if re.search('.csv$', f):
                            with open(mypath+'/'+x+'/'+f, 'r') as f:
                                reader = csv.reader(f, delimiter=',')
                                next(reader, None) #skip headers
                                #headers = csv.DictReader(f, delimiter=',')
                                for line in reader:
                                    #for csv
                                    csv_writer = csv.writer(csvfile, delimiter=',')  
                                    csv_writer = csv_writer.writerow(line)  
                                    longitude.append(line[1]) 
                                    latitude.append(line[0])
                                    thickness.append(line[3])
                                    bedrock.append(line[7])
                                    surface.append(line[6])
                                
    #save csvfile (for QGis for example) 
    print('Saving CSV file....')
                                      
    #change list in array      
    print('changing list to numpy array...')
    longitude = np.asarray(map(float, longitude[0:np.size(longitude)]))
    print('    longitude ok')
    latitude = np.asarray(map(float, latitude[0:np.size(latitude)]))
    print('    latitude ok' )
    thickness = np.asarray(map(float,thickness[0:np.size(thickness)]))
    print('    thickness ok')
    bedrock = np.asarray(map(float,bedrock[0:np.size(bedrock)]))
    print('    bedrock ok')
    surface = np.asarray(map(float,surface[0:np.size(surface)]))
    print('    surface ok')
    
    print('Lecture done. Saving Data in .npy format...')
    MATRIX = np.array([longitude, latitude, thickness, bedrock, surface])
    MATRIX = np.transpose(MATRIX)
    #np.savetxt('matrix.txt', MATRIX, delimiter=' ', newline='\n' ) 
    np.save('data_'+filename, MATRIX)
               
    #return longitude, latitude, thickness, bedrock, surface                                       

def lecture_airborn_pre_ice_bridge(*args):#(mypath1, mypath2):
    longitude = []
    latitude = []
    thickness = []
    bedrock = []
    surface = []
    #folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 

    with open('Pre_Ice_Bridge_NASA.csv', 'w') as csvfile: 
        for mypath in args:
            print('Main Directory:', mypath)
            folders = [ x for x in listdir(mypath) if isdir(join(mypath,x)) ] 
            for x in folders:
                if re.search('^2002', x) or re.search('^2004', x):# or re.search('^2011', x) :
                    print('    start reading directory: ', x)
                    onlyfiles = [ f for f in listdir(mypath+'/'+x) if isfile(join(mypath+'/'+x,f)) ] 
                    for f in onlyfiles :
                        if re.search('.csv$', f):
                            with open(mypath+'/'+x+'/'+f, 'r') as f:
                                reader = csv.reader(f, delimiter=',')
                                next(reader, None) #skip headers
                                #headers = csv.DictReader(f, delimiter=',')
                                for line in reader:
                                    #for csv
                                    csv_writer = csv.writer(csvfile, delimiter=',')  
                                    csv_writer = csv_writer.writerow(line)  
                                    longitude.append(line[1]) 
                                    latitude.append(line[0])
                                    thickness.append(line[4])
                                    bedrock.append(line[7])
                                    surface.append(line[6])
                                
    #save csvfile (for QGis for example) 
    print ('Saving CSV file....')
                                      
    #change list in array      
    print('changing list to numpy array...')
    longitude = np.asarray(map(float, longitude[0:np.size(longitude)]))
    print('    longitude ok')
    latitude = np.asarray(map(float, latitude[0:np.size(latitude)]))
    print('    latitude ok' )
    thickness = np.asarray(map(float,thickness[0:np.size(thickness)]))
    print('    thickness ok')
    bedrock = np.asarray(map(float,bedrock[0:np.size(bedrock)]))
    print('    bedrock ok')
    surface = np.asarray(map(float,surface[0:np.size(surface)]))
    print('    surface ok')
    
    print('Lecture done. Saving Data in .npy format...')
    MATRIX = np.array([longitude, latitude, thickness, bedrock, surface])
    MATRIX = np.transpose(MATRIX)
    #np.savetxt('matrix.txt', MATRIX, delimiter=' ', newline='\n' ) 
    np.save('Pre_Ice_Bridge_NASA', MATRIX)
               
    return longitude, latitude, thickness, bedrock, surface                               
                                                
                                                                        
                                                                                                                        
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km
    
    
# Data manipulation:

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments    
    
def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc    
 
            
            
                      
            