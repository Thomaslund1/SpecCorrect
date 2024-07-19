"""Global Packadges"""
import numpy as np
import h5py as hpy
import pandas as pd 


def interpolate(wls):
    """
    Function to interpolate nans with +-10 point median
    @param wls : array/list
        >1d list or array of wavelengths with nans that need interpolating
    @returns wls : 
        same array or list but with interpolated nans
    """
    mask = np.isnan(wls)

    # Find the indices of NaN values
    nan_indices = np.where(mask)

    # Replace NaN values with the median of surrounding 10 points
    for i, j in zip(*nan_indices):
        if j == 0:
            wls[i, j] = 5000
        else:
            wls[i, j] = np.nanmedian(wls[i][j-10:j+10])
    return wls


def abs2OrdInd(absPos,order,index):
    """
    a function to return the order and relative index from an absolute index value
    @param absPos : int
        The absolute position that we want to turn to an orderred pair
    @param order : list/array
        the list of orders for the same measurment as the absPos value
    @param index : list/array
        the list of indexes for the same measurment as the abPos value
    """
    return([order[absPos],index[absPos])

def order_index_selector (orders, indices, order, index, ref):
    """
    a function which takes in a selected order, index value, and reference file
    @ param orders : array 
    takes the desired orders array from a HDF5 file 
    @param indices : array 
    takes the desired indices array from a HDF5 file 
    @param order: int 
    takes the desired order number 
    @param index : int 
    takes the desired index number
    @param ref : int
    takes the integer value of the desired reference file, should not matter since the indices and orders should be the same across file choices. 
    @returns 
    a list of absolute index values (absolute index being the index within the hdf5 file structure from 0 - 3575 ish).
    """
    abs_indices = [] 
    for i in range(len(orders[ref])):
        if orders[ref][i] == order and indices[ref][i] == index:
            abs_indices.append(i)       
    return abs_indices

"""Index Binning"""
def group_by_num(num, data_loc, ref, wl, pix, order_indices, pix_array = None, indices_array = None):        
        """
    Groups data by index into sequnetial bins. Returns index values for each bin and if desired lined up wavelength and pixel data. 
    @param num : int
        number of entries per bin
    @param data_loc : str / array  
        file path for data to chunk. If one wants to insert a specific array instead, put the array name in for data_loc instead.
        Note: If arrays are used instead of filepath, pix_arrays and indices_arrays must also be filled with the their respective 
        array names. 
    @param ref : int
        which measurment to use as the reference date for pulling indexes
    @param wl : bol 
        Set to 1 if outputted wavelengths in specific bins is desired. 
    @param pix : bol
        Set to 1 if outputted pixels in specific bins is desired.
    @param order_indecies : bol
       Set to 1 if outputted indices in specific bins is desired.
    @param pix_array : array 
        Default is fot this argument to be ignored. If one wants to use already loaded arrays instead of pointing to a file dircetory 
        an array of pixel data must be pointed to. 
    @param indices_array : array 
        Default is for this argument to be ignored.  If one wants to use already loaded arrays instead of pointing to a file dircetory 
        an array of indices must be pointed to.
    @returns list_ind, list_wavl, list_pix, list_ord_ind
        Lists of the indices/wavleneght/pix/order binned by the desired number. 
    """

        if type(data_loc) == str:
            cenM_ar = np.array(hpy.File(data_loc + 'All_centroidWl.hdf5', 'r')['dat'][:]) # points at specific location 
            cenM_pix = np.array(hpy.File(data_loc + 'All_centroidPix.hdf5', 'r')['dat'][:]) # points at specific cenM_pix
            indices = np.array(hpy.File(data_loc + 'All_index.hdf5', 'r')['dat'][:]) # points at specific indices
        else:
            cenM_ar = data_loc
            cenM_pix = pix_array
            indices = indices_array

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
            list_pix = [] 
            for i in range(len(cenM_ar[ref])):
                if i%num == 0 and i != 0:
                    sublist_pix = [] 
                    for i in range(i-num, i):
                        sublist_pix.append(cenM_pix[ref][i])
                    list_pix.append(sublist_pix)

        if order_indices == 1: 
            list_ord_ind = [] 
            for i in range(len(indices[ref])):
                if i%num == 0 and i != 0:
                    sublist_ordind = [] 
                    for i in range(i-num, i):
                        sublist_ordind.append(indices[ref][i])
                    list_ord_ind.append(sublist_ordind)
    
        return list_ind, list_wavl, list_pix, list_ord_ind


"""order Stuff"""
def getOrder(file, num, bins, dirFile):
    """
    Helper funtion to groupByOrder to chunk one order
    @param file : array/list
        the full list of orders as pulled from the hdf5
    @param num : int
        which order to analyze
    @param bins : int
        how many bins to chop into
    @param dirFile : 
        a data file to split by the index values calculated
    @returns out : np.array
        a 2d array with an entry for each bin that contains all the data chuncked into it
    """
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

def groupByOrder(orders, data, bins):
    """
    Proocesses the set of given data by order bins
    @param orders : list/array
        the list of orders from a hdf5 file
    @param data : list/array
        the data to chunk as pulled from hdf5 file
    @param bins : int
        how many bins you want returned
    @returns listOfBins : list
        the chopped data in the format [(orders),(bins),(measurments by date),(pixel/index)]
        i.e listOfBins[15,2,5,6] returns the 16th order, 3rd bin, 6th measruement, 7th index value of whatever was passed by data
    """
    listOfBins = []
    
    for q in range(abs(int(orders[0][0]) - int(orders[0][-1])) - 1): # calculates first minus last order number to know how many orders to process
        print(f"Processing order {int(q+1)}/{(abs(int(orders[0][0]) - int(orders[0][-1])) - 1)}", end="\r") # prints an estimated progress report based on the above calculation
        order_number = q + orders[0][0] + 1 # adjusts for order offsets like NEID Etalon starting at order 26
        
        dat = getOrder(orders, order_number, bins, data)
        newLis = []
        
        for i in range(len(dat)): # reformats the binned measurments and adds them to a constant list
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

"""Velocity Function"""
def wl2vel(wls,ref=50):
    """
    converts wavelength values to velocities by a reference date
    @param wls : list/array
        the wavelengths to convert
    @param ref : int (default 50)
        the reference date to convert to velocities from (getting change in wl)
    @returns out : list
        the list of velocities
    """
    reference_row = wls[ref].copy()
    out = []
    for i in range(len(wls)):
        out.append((np.subtract(wls[i],reference_row)/reference_row)*299792458)
    return(out)

def getVels(wavelengths,orders,bins,ordVsInd):
    """
    repackadged version of both functions to return a binned list of velocities
    @param wavelengths : list/array
        the list of wavelengths as pulled from the hdf5 file
    @param orders : list/array
        list of orders as pulled from hdf5
    @param bins : int
        how many bins per order if using per order, or how many index values if using by index
    @param ordVsInd : bol
        0 = bin by order | 1 = bin by index
    @returns out : list
        the sliced data in terms of velocities 
    """
    vels = wl2vel(wavelengths)
    if(not ordVsInd):
        return groupByOrder(orders,vels,bins)
    else:
        return group_by_num(bins, vels, 50, 1)
        

def filter_vels(data, indices, ref):
    """"
    function used for binning velocities by desired indices
    @param data : array 
        an array of velocities 
    @param : array/list
        an array or list of desired indices 
    @param : int  
        an integer which points to a specific reference file to use 
    @returns out: list 
        velocities filtered by the desired indices.
    """""

    out = []
    data1 = np.array(data,dtype = np.float64)
    indices = np.array(indices,dtype = "int")
    for i in range(len(indices)):
        out.append(data1[ref][indices[i]])
    return out



def abs_ind_vels(data, abs_ind):
    """
    function used for calculating the velocities of a specific absolute index value calculated from desired index and order 
    @param : array 
        an array of wavlenghts (or pixles) that one wants to parse through to calculate velocities. 
    @param : abs_ind 
        an integer value that points the code to a desired absolute index (0 - 3575)
    @returns out : list 
        velocities for just the specified absolute index which should from above correspond with a (order, index) ordered pair. 
    """
    c = 299792458
    velocities = [] 
    for i in range(len(data)):
        dif = data[i][abs_ind] - data[0][abs_ind]
        div = dif/data[0][abs_ind]
        vel = div/c
        velocities.append(vel)
    return velocities



"""CSV Function"""
def out_csv(file_path,out):
    """
    Outputs velocities data into a csv file. 
    @param: file_path : str 
        Note: the file path cannot just be a directory it should include the desired file directory with the intended name for the file at the end of 
        the directory. 
    @param out : list/array
        data to store in csv
    
    """
   
    data = {
        "Velocities": out
    }

    df = pd.DataFrame(data)
    df.to_csv(file_path, index=False)

    return
