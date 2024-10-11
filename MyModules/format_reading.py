 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Created on Wed Jan 16 16:54:52 2019

@author: cmosbeux

Description :   vtu     :  Read a vtu parallel file and give desired output 
                           if they exist.
                netcdf  :  Read netcdf and return vectors x,y and variable
                
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

import re
from os import listdir
from os.path import isfile, join
import vtk
from vtk.util import numpy_support as VN
import numpy as np
from netCDF4 import Dataset


def vtu(input_directory, file_name, list_of_var):
    '''
    input: 
        -input_directory : directory containing files 
        -file_name : generic name of the file (part of it allowing identification)
        -list_of_var: list of var wanted in the output "dictionnary list"
    
    output:
        -coord_x
        -coord_y
        -dictionnary containing all the variables, e.g. v['your_variable']
    '''
    
    onlyfiles = [ f for f  in listdir(input_directory) 
        if isfile(join(input_directory+'/'+f)) ]

    
    onlyfiles=sorted(onlyfiles)
    
    #What variables do you want to plot/use?
    coords = []
    v = {}
    for i in list_of_var:
        v[i]=[]
    
    try:
        k=0
        for f in onlyfiles:
            #check for "root" name 
            fname = re.sub('\.vtu$', '', f)
            input_name = re.sub('\.vtu$', '', file_name)
            fext= f.split('.')[-1]
            if (fext!='vtu'):
                continue
            flag = True
            #check for partitions
            for fpart in input_name.split('*'):
                #second condition fpart !='' because of split in py3
                if fpart != '':
                    if fpart in fname:
                        continue
                    else:
                        flag = False
            if not flag:
                continue
            #Read the source file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(input_directory+'/'+f)
            reader.Update()
            output=reader.GetOutput()
            PointData=output.GetPointData()
            coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))

            #concatenate all the partitions
            for i in list_of_var:
                array=VN.vtk_to_numpy(PointData.GetArray(i))
                v[i].extend(array)
            k+=1
            
    except AttributeError:
        print('variable', i, 'does not exist...')
        
    except IndexError:
        print('Problem with the file (it seems empty)... check path and name.')
        
    print('VTU (%d partitions) %s inspected...' % (k, file_name))
            
    #change lists to array
    def arr(i):
        i=np.asarray(i)
        return i
    #coords
    coord_x=arr(coords).T[0]
    coord_y=arr(coords).T[1]
    
    return coord_x, coord_y, v


def vtu_cells(input_directory, file_name, list_of_var, data_type='point'):
    '''
    input:
        - input_directory: directory containing files
        - file_name: generic name of the file (part of it allowing identification)
        - list_of_var: list of variables wanted in the output "dictionary list"
        - data_type: type of data to extract ('point' or 'cell')
    
    output:
        - coord_x
        - coord_y
        - dictionary containing all the variables, e.g. v['your_variable']
    '''
    
    onlyfiles = [f for f in listdir(input_directory) if isfile(join(input_directory, f))]
    onlyfiles = sorted(onlyfiles)

    # What variables do you want to plot/use?
    coords = []
    v = {}
    for i in list_of_var:
        v[i] = []
    
    try:
        k = 0
        for f in onlyfiles:
            # check for "root" name 
            fname = re.sub('\.vtu$', '', f)
            input_name = re.sub('\.vtu$', '', file_name)
            fext = f.split('.')[-1]
            if fext != 'vtu':
                continue
            flag = True
            # check for partitions
            for fpart in input_name.split('*'):
                if fpart != '':
                    if fpart in fname:
                        continue
                    else:
                        flag = False
            if not flag:
                continue
            
            # Read the source file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(join(input_directory, f))
            reader.Update()
            output = reader.GetOutput()
            
            if data_type == 'point':
                PointData = output.GetPointData()
                coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))
            elif data_type == 'cell':
                CellData = output.GetCellData()
                num_cells = output.GetNumberOfCells()
                coords.extend([(output.GetCell(i).GetPoints().GetPoint(0)[0], output.GetCell(i).GetPoints().GetPoint(0)[1]) for i in range(num_cells)])
            
            # Concatenate all the partitions
            for i in list_of_var:
                if data_type == 'point':
                    array = VN.vtk_to_numpy(PointData.GetArray(i))
                elif data_type == 'cell':
                    array = VN.vtk_to_numpy(CellData.GetArray(i))
                v[i].extend(array)
            k += 1
            
    except AttributeError:
        print('Variable', i, 'does not exist...')
        
    except IndexError:
        print('Problem with the file (it seems empty)... check path and name.')
        
    print('VTU (%d partitions) %s inspected...' % (k, file_name))
            
    # Change lists to array
    def arr(i):
        i = np.asarray(i)
        return i
    

    # Coords
    if data_type == 'point':
        coord_x = arr(coords).T[0]
        coord_y = arr(coords).T[1]
    else:
       coord_x = [c[0] for c in coords]
       coord_y = [c[1] for c in coords]
        
    
    return coord_x, coord_y, v


def vtu3D(input_directory, file_name, list_of_var):
    '''
    input: 
        -input_directory : directory containing files 
        -file_name : generic name of the file (part of it allowing identification)
        -list_of_var: list of var wanted in the output "dictionnary list"
    
    output:
        -coord_x
        -coord_y
        -dictionnary containing all the variables, e.g. v['your_variable']
    '''
    
    onlyfiles = [ f for f  in listdir(input_directory) 
        if isfile(join(input_directory+'/'+f)) ]

    
    onlyfiles=sorted(onlyfiles)
    
    #What variables do you want to plot/use?
    coords = []
    v = {}
    for i in list_of_var:
        v[i]=[]
    
    try:
        k=0
        for f in onlyfiles:
            fname = re.sub('\.vtu$', '', f)
            input_name = re.sub('\.vtu$', '', file_name)
            fext= f.split('.')[-1]
            if (fext!='vtu'):
                continue
            flag=True
            for fpart in input_name.split('*'):
                if fpart in fname :
                    continue
                else:
                    flag=False
            if not flag:
                continue
            
            #Read the source file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(input_directory+'/'+f)
            reader.Update()
            output=reader.GetOutput()
            PointData=output.GetPointData()
            coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))
            k+=1
            
            for i in list_of_var:
                    array=VN.vtk_to_numpy(PointData.GetArray(i))
                    v[i].extend(array)
                    
    except IndexError:
        print('Problem with the file (it seems empty)... check path and name.')   
        
    except AttributeError:
        print('variable', i, 'does not exist...')
                
        
    print(('VTU (%d partitions) %s inspected...') % (k, file_name))
    print(('\t dimension = (%d nodes ,%d axis)') % np.array(coords).shape)
    #change lists to array
    def arr(i):
        i=np.asarray(i)
        return i
    #coords
    coord_x=arr(coords).T[0]
    coord_y=arr(coords).T[1]
    coord_z=arr(coords).T[2]
    return coord_x, coord_y, coord_z, v
    

def vtu3D_bis(input_directory, file_name, list_of_var, var_type='point'):
    '''
    input: 
        -input_directory : directory containing files 
        -file_name : generic name of the file (part of it allowing identification)
        -list_of_var: list of var wanted in the output "dictionary list"
        -var_type: point (default) or cell 
    
    output:
        -coord_x
        -coord_y
        -coord_z
        -dictionnary containing all the variables, e.g. v['your_variable']
        -face_connections: list of face connections for triangulation
    '''
    onlyfiles = [f for f in listdir(input_directory) if isfile(join(input_directory, f))]
    onlyfiles = sorted(onlyfiles)
    
    coords = []
    coords_cell = []
    v = {}
    for i in list_of_var:
        v[i] = []

    
    face_connections = []
    
    try:
        k = 0
        for f in onlyfiles:
            fname = re.sub('\.vtu$', '', f)
            input_name = re.sub('\.vtu$', '', file_name)
            fext = f.split('.')[-1]
            if fext != 'vtu':
                continue
            flag = True
            for fpart in input_name.split('*'):
                if fpart in fname:
                    continue
                else:
                    flag = False
            if not flag:
                continue
            
            # Read the source file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(join(input_directory, f))
            reader.Update()
            output = reader.GetOutput()
            
            #get point data
            PointData = output.GetPointData()
            coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))
            cells = output.GetCells()
            cell_array = VN.vtk_to_numpy(cells.GetData())
            num_cells = cells.GetNumberOfCells()
            cell_sizes = cell_array[:num_cells]
            cell_offsets = np.cumsum(cell_sizes)
            cell_connectivity = cell_array[cell_offsets[-1]:]
            face_connections.extend(cell_connectivity)
            
            # Concatenate all the partitions
            if var_type == 'point':
                for i in list_of_var:
                    array = VN.vtk_to_numpy(PointData.GetArray(i))
                    v[i].extend(array)
              
            #get cell data
            if var_type == 'cell':
                CellData = output.GetCellData()
                coords_cell.extend([(output.GetCell(i).GetPoints().GetPoint(0)[0], output.GetCell(i).GetPoints().GetPoint(0)[1], output.GetCell(i).GetPoints().GetPoint(0)[2]) for i in range(num_cells)]) 
        
                # Concatenate all the partitions
                for i in list_of_var:
                    array = VN.vtk_to_numpy(CellData.GetArray(i))
                    v[i].extend(array)
                
            k += 1
            
        if var_type == 'cell':
            coords = coords_cell        
    
    except IndexError:
        print('Problem with the file (it seems empty)... check path and name.')   
        
    except AttributeError:
        print('Variable', i, 'does not exist...')
                
    print(('VTU (%d partitions) %s inspected...') % (k, file_name))
    print(('\t dimension = (%d nodes, %d axis)') % np.array(coords).shape)
    
    # Change lists to array
    def arr(i):
        i = np.asarray(i)
        return i
    
    # Extract coordinates
    coord_x = arr(coords).T[0]
    coord_y = arr(coords).T[1]
    coord_z = arr(coords).T[2]
    
 
    return coord_x, coord_y, coord_z, v, face_connections
    


def vtu_boundaries(input_directory, file_name, list_of_var):
    '''
    input: 
        -input_directory : directory containing files 
        -file_name : generic name of the file (part of it allowing identification)
        -list_of_var: list of var wanted in the output "dictionnary list"
    
    output:
        -coord_x
        -coord_y
        -dictionnary containing all the variables, e.g. v['your_variable']
    '''
    
    onlyfiles = [ f for f  in listdir(input_directory) 
        if isfile(join(input_directory+'/'+f)) ]

    
    onlyfiles=sorted(onlyfiles)
    
    #What variables do you want to plot/use?
    coords = []
    v = {}
    for i in list_of_var:
        v[i]=[]
    
    try:
        k=0
        for f in onlyfiles:
            #check for "root" name 
            fname = re.sub('\.vtu$', '', f)
            input_name = re.sub('\.vtu$', '', file_name)
            fext= f.split('.')[-1]
            if (fext!='vtu'):
                continue
            flag = False
            #check for partitions
            for fpart in input_name.split('*'):
                #second condition fpart !='' because of split in py3
                if fpart in fname and fpart != '':
                    flag = True
                    continue
                else:
                    continue
            if not flag:
                continue
            #Read the source file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(input_directory+'/'+f)
            reader.Update()
            output=reader.GetOutput()
            PointData=output.GetPointData()
            CellData=output.GetCellData()
            coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))

            #concatenate all the partitions
            for i in list_of_var:
                array=VN.vtk_to_numpy(PointData.GetArray(i))
                v[i].extend(array)
            k+=1
            
    except AttributeError:
        print('variable', i, 'does not exist...')
        
    except IndexError:
        print('Problem with the file (it seems empty)... check path and name.')
        
    print('VTU (%d partitions) %s inspected...' % (k, file_name))
            
    #change lists to array
    def arr(i):
        i=np.asarray(i)
        return i
    #coords
    coord_x=arr(coords).T[0]
    coord_y=arr(coords).T[1]
    
    return coord_x, coord_y, v


#%%read NetCDF
class netcdf:
    
    @staticmethod
    def readline(fname,varname, xx='x', yy='y', format='NETCDF3_CLASSIC'):
        '''
        input: 
            -fname  : file name (+path) 
            -varname: variable to load
            -format : Default is 'NETCDF3_CLASSIC'
        output:
            -coord_x 1D array
            -coord_y 1D array
            -variable 1D array
        '''
        dataset1=Dataset(fname, 'r',  format=format) 
        x = dataset1.variables[xx][:]
        y = dataset1.variables[yy][:]
        v = dataset1.variables[varname][:]
        v= v.flatten()
        xx,yy,vv =[],[],[]
        k=0
        for i in x:
            for j in y:        
                xx.append(j)
                yy.append(i) 
                vv.append(v[k])
                k+=1
        return xx,yy,vv
    
    @staticmethod
    def readgrid(fname,varname, xx='x', yy='y', format='NETCDF3_CLASSIC'):
        '''
        input: 
            -fname : file name (+path) 
            -varname : variable to load
            -format : formating
        output:
            -coord_x 1D array
            -coord_y 1D array
            -variable 2D array
        '''
        dataset1=Dataset(fname, 'r',  format=format) 
        x = dataset1.variables[xx][:]
        y = dataset1.variables[yy][:]
        v = dataset1.variables[varname][:]
        return x,y,v
    
    @staticmethod
    def write(fname,vardic, xname='x',yname='y', unit='-', description = None, format='NETCDF3_CLASSIC'):
        '''write a netcdf with the diferent variables from a collection/dic
        containing x (1D), y (1D) and variables (2D)'''
        dataset=Dataset(fname, 'w',  format=format) 
        m,n=len(vardic[xname]),len(vardic[yname])
        x=dataset.createDimension(xname,m)
        y=dataset.createDimension(yname,n)
        coord_y = dataset.createVariable(yname,"f4",(yname,))
        coord_x = dataset.createVariable(xname,"f4",(xname,))
        coord_x[:]=vardic[xname]
        coord_y[:]=vardic[yname]
        k=0
        varsave={}
        for i in vardic:
            if i!=xname and i!=yname:
                varsave[i]=dataset.createVariable(i,"f4",(yname,xname,), fill_value=-1e8)
                varsave[i][:,:]=vardic[i]
                varsave[i].units = unit
                k+=1
        coord_x.units = 'm'
        coord_y.units = 'm'
        
        if description is not None:
            dataset.description = description
        dataset.close()
           
        
