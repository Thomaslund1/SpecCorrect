
"""Global Packadges"""
import numpy as np
from datetime import datetime
import h5py as hpy
import sys
import matplotlib.dates as mdates
import pandas as pd

# function for masking out the nan values from the data set.
def interpolate(wls):
    mask = np.isnan(wls)

    # Find the indices of NaN values
    nan_indices = np.where(mask)

    # Replace NaN values with the previous non-NaN value along each row
    for i, j in zip(*nan_indices):
        if j == 0:
            wls[i, j] = 5000
        else:
            wls[i, j] = np.nanmedian(wls[i][j-10:j+10])
    return wls

"""Index Binning"""
#for the function to work, data_loc must be the directory where the consolidated data is stored. For the arguments wl, pix, and order_indices 
#enter 1 if you want those lists returned. 
def group_by_num(num, data_loc, ref, wl, pix, order_indices):    
    
    cenM_ar = np.array(hpy.File(data_loc + 'All_centroidWl.hdf5', 'r')['dat'][:]) # points at specific location 

    list_ind = [] 

    for i in range(len(cenM_ar[ref])):
        if i%num == 0 and i != 0:
            sublist_ind = [] 
            for i in range(i-num, i):
                sublist_ind.append(i)
            list_ind.append(sublist_ind)
    
    if wl == 1:
        list_wavl = [] 
        for i in range(len(cenM_ar[ref])):
            if i%num == 0 and i != 0:
                sublist_wavl = [] 
                for i in range(i-num, i):
                    sublist_wavl.append(cenM_ar[ref][i])
                list_wavl.append(sublist_wavl)

    if pix == 1:
        cenM_pix = np.array(hpy.File(data_loc + 'All_centroidPix.hdf5', 'r')['dat'][:])
        list_pix = [] 
        for i in range(len(cenM_ar[ref])):
            if i%num == 0 and i != 0:
                sublist_pix = [] 
                for i in range(i-num, i):
                    sublist_pix.append(cenM_pix[ref][i])
                list_pix.append(sublist_pix)

    if order_indices == 1: 
        list_ord_ind = [] 
        indices = np.array(hpy.File(data_loc + 'All_index.hdf5', 'r')['dat'][:])
        for i in range(len(indices[ref])):
            if i%num == 0 and i != 0:
                sublist_ordind = [] 
                for i in range(i-num, i):
                    sublist_ordind.append(indices[ref][i])
                list_ord_ind.append(sublist_ordind)
    
    return list_ind, list_wavl, list_pix, list_ord_ind

"""order Stuff"""
def getOrder(file, num, bins, dirFile):
    mask = (file == num)
    temp = []
    for i in range(len(dirFile)):
        temp.append(dirFile[i][mask[i]])
    filtered_data = np.array(temp) 
    sizeOfBin = filtered_data.shape[1] // bins  # Size of each full bin
    
    trimmed_data_length = sizeOfBin * bins  # Length to trim to
    trimmed_data = filtered_data[:, :trimmed_data_length] 
    
    out = np.split(trimmed_data, bins, axis=1)
    return np.array(out)

def groupByOrder(orders, wls, bins):
    eq_rad_velocities = wl2vel(wls)
    listOfBins = []
    

    for q in range(abs(int(orders[0][0]) - int(orders[0][-1])) - 1):
        print(f"Processing order {int(q+1)}/{(abs(int(orders[0][0]) - int(orders[0][-1])) - 1)}", end="\r")
        order_number = q + orders[0][0] + 1
        
        dat = getOrder(orders, order_number, bins, eq_rad_velocities)
        newLis = []
        
        for i in range(len(dat)):
            subl = []
            for j in range(len(dat[i])):
                subl.append(dat[i][j])
            newLis.append(subl)
        newList = np.array(newLis)
        
        listOfEras = []
        for i in range(bins):
                listOfEras.append(newList[i])
        listOfBins.append(listOfEras)
    
    return listOfBins

"""Velocity Functions"""
def wl2vel(wls):
    reference_row = wls[0].copy()
    return(((wls-reference_row)/reference_row)*299792458)

#calculates velocities from a selected order with the specified number of bins. 
def getVels(wavelengths,orders,bins):
    waves = wl2vel(wavelengths)
    return groupByOrder(orders,waves,bins)


#given vels (or any data of your choice), this function filters through the velocities based on provided indices ranges. A reference file is also required as it is the file from which the velocities will be pulled from. 
def filter_vels(data, indices, ref):
    out = []
    data1 = np.array(data,dtype = np.float64)
    indices = np.array(indices,dtype = "int")
    for i in range(len(indices)):
        out.append(data1[ref][indices[i]])
    return out
