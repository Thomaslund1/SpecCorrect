"""Global Packadges"""
import numpy as np
import h5py as hpy
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.signal import lombscargle
import astropy
import copy

"""Reference File Generator"""
# note the function input must be an npy file
def reference_file_maker(file_path,name):
    """
    Function used to generate new reference files to run linefits data through. 
    @param file_path : str 
    The path of the npy file being used to generate a new reference file. 
    @param name : str 
    name of the new reference file that will be generated.
    @return file
    A new saved npy reference file to use on additional linefits runs. 
    """
    #loading in the file that will be used to generate a reference file
    file = np.load(file_path, allow_pickle = 1)[()]
    #creating a copy of the file and storing it in the variable fi 
    fi = copy.deepcopy(file)

    """extracting the data corresponding with the 'centroid_pix' key. This
    is the only data point used when generating a new reference files."""

    for i in file.keys():
        for j in file[i].keys():
            fi[i][j] = file[i][j]['centroid_pix']

    np.save(name,fi)
    return 

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
            wls[i, j] = np.nanmedian(wls[i][j - 10 : j + 10])
    return wls


def boxcar_median(data, window_size):
    """
    takes the moving median of a 1d set of data
    @param data : list/array
        the data to smooth
    @param window_size : int
        how many points per median should be included
    """
    half_window = window_size // 2
    smoothed_data = np.zeros(len(data), dtype=np.float64)

    for i in range(len(data)):
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        smoothed_data[i] = np.median(data[start:end])

    return smoothed_data


"""
It is about time
"""


def convert_time(data_list):
    """
    converts a list of times into a list of datetime values
    @param data list : list/array
    the time data that one wants to convert to datetime values
    @returns time : list
    a list of datetime values generated from the list of times from data_list
    """

    time2 = (
        []
    )  # creating a list to store all the times in a way that time date can read

    for i in data_list:
        a_str = "%i" % i
        a_str = a_str[:8] + "T" + a_str[8:]
        time2.append(a_str)

    time = (
        []
    )  # creating a list that stores the datetime data. Serves as x-axis on plots

    format_data = "%Y%m%dT%H%M%S"  # Format string matching the time data structure

    for i in range(len(time2)):
        date = datetime.strptime(time2[i], format_data)
        time.append(date)

    return time


#binning signifies if one wants time in seconds, minutes, hours, days, or weeks. seconds = 0, minutes = 1, ... 
def datetime_conversionify(time_list, binning):
    datetime_objects_seconds = [] 

    if binning == 0:
        for i in range(len(time_list)):
            a_timedelta = time_list[i] - time_list[0]
            seconds = a_timedelta.total_seconds()
            datetime_objects_seconds.append(seconds)
            
    elif binning == 1:
        seconds_to_minutes = []
        for i in range(len(time_list)):
            a_timedelta = time_list[i] - time_list[0]
            seconds = a_timedelta.total_seconds()
            seconds_to_minutes.append(seconds)
        for i in range(len(seconds_to_minutes)):
            minute = (seconds_to_minutes[i]/60)
            datetime_objects_seconds.append(minute)

    elif binning == 2: 
        seconds_to_minutes = []
        minutes_to_hours = [] 
        for i in range(len(time_list)):
            a_timedelta = time_list[i] - time_list[0]
            seconds = a_timedelta.total_seconds()
            seconds_to_minutes.append(seconds)
        for i in range(len(seconds_to_minutes)):
            minute = (seconds_to_minutes[i]/60)
            minutes_to_hours.append(minute)
        for i in range(len(minutes_to_hours)):
            hour = (minutes_to_hours[i]/60)
            datetime_objects_seconds.append(hour)

    elif binning == 3: 
        seconds_to_minutes = []
        minutes_to_hours = [] 
        hours_to_days = [] 
        for i in range(len(time_list)):
            a_timedelta = time_list[i] - time_list[0]
            seconds = a_timedelta.total_seconds()
            seconds_to_minutes.append(seconds)
        for i in range(len(seconds_to_minutes)):
            minute = (seconds_to_minutes[i]/60)
            minutes_to_hours.append(minute)
        for i in range(len(minutes_to_hours)):
            hour = (minutes_to_hours[i]/60)
            hours_to_days.append(hour)
        for i in range(len(hours_to_days)):
            day = (hours_to_days[i]/24)
            datetime_objects_seconds.append(day)

    elif binning == 4: 
        seconds_to_minutes = []
        minutes_to_hours = [] 
        hours_to_days = [] 
        days_to_weeks = [] 
        for i in range(len(time_list)):
            a_timedelta = time_list[i] - time_list[0]
            seconds = a_timedelta.total_seconds()
            seconds_to_minutes.append(seconds)
        for i in range(len(seconds_to_minutes)):
            minute = (seconds_to_minutes[i]/60)
            minutes_to_hours.append(minute)
        for i in range(len(minutes_to_hours)):
            hour = (minutes_to_hours[i]/60)
            hours_to_days.append(hour)
        for i in range(len(hours_to_days)):
            day = (hours_to_days[i]/24)
            days_to_weeks.append(day)
        for i in range(len(days_to_weeks)):
            week = (days_to_weeks[i]/7)
            datetime_objects_seconds.append(week)
            

    return datetime_objects_seconds


def CompatibleDataArraysIND(time_data_small, time_data_large):
    """
    Function intended to be used on HPF data where the science fiber data and calibration fiber do not lineup.
    @param time_data_small
    the array of time data from the hdf5 file. This array is meant to be the data from the fiber source with
    less time data.
    @param time_data_large
    the array of data from the hdf5 file. This array is meant to be the data from the fiber source with more
    time data.
    @returns compat_indices
    A list of indices where the recorded times of the measurments from each file tend to lineup better.
    """
    # converting the time arrays to datetimes using the convert_time function from above.
    small_time = convert_time(time_data_small)
    large_time = convert_time(time_data_large)

    # named minimum_indices because it will keep track of the index values where the difference between the
    # recorded cal fiber time and sci fiber time are at a minimum

    compat_indices = []
    for j in range(len(small_time)):
        magnitudes = []
        subtract = []
        for i in range(len(large_time)):
            # finding the difference in seconds between time points in the larger array and smaller array
            diff = large_time[i] - small_time[j]
            subtract.append(diff.total_seconds())
        for i in range(len(subtract)):
            # taking the absolute value of the differences so that the true minimum magnitude can be found
            absolute = abs(subtract[i])
            magnitudes.append(absolute)

        # finding the minimum of the magnitudes, where the minimum index is, and putting it in the compat_indices array
        minimum = np.min(magnitudes)
        min_ind = magnitudes.index(minimum)
        compat_indices.append(min_ind)
        # progress bar
        print(str(j) + "/" +str(len(time_data_small)), end="\r")
    return compat_indices


def CompatibleWavelengths(wavl_data, compat_indices):
    """
    Function intended to create 'lined up' values for centroid data across scical fibers that potentially have different amounts of data.
    @param wavl_data
    the function is intended to take data from a larger array than the array the data will be compared to. For example,
    for sci-cal differentials, if cal has much more data that is misaligned with sci, cal centroid data should be passed in the
    argument for wavl_data.
    @param compat_indices
    the indices generated from the CompatibleArraysIND.
    @returns compatible_cen
    A list of centroid wavelenghts that should be compatible with the smaller arrays
    """
    compatible_cen = []
    # adding the compatible wavelength centroid data as determined by the CompatibleDataArraysIND function
    for i in compat_indices:
        compatible_cen.append(wavl_data[i])
    return compatible_cen


def lSPerodogram(periods, time, flux, graph=0):
    """
    takes the lomb scargle periodogram of the given data to find periodicities
    @param periods : list/array
        a set of ints/floats that the function will try to find periods of in the data, recomeded to use np.logscale
    @param time : list/array
        the time axis of data to analyze
    @param flux : list/array
        the intensities that correspond to the time data
    @param graph : bol
        whether to graph the normalized periodogram of the base data
    @return power : list
        the relative contributions of each input period to the given data
    """

    frequencies = 1 / periods
    power = lombscargle(time, flux, frequencies)

    if graph:
        plt.figure(figsize=(10, 4))
        plt.plot(periods, power)
        plt.title("Lomb-Scargle Periodogram")
        plt.xlabel("Period")
        plt.ylabel("Power")
        plt.xscale("log")
        plt.grid(True)
        plt.gca().invert_xaxis()
        plt.show()
    return power


def abs2OrdInd(absPos, order, index):
    """
    a function to return the order and relative index from an absolute index value
    @param absPos : int
        The absolute position that we want to turn to an orderred pair
    @param order : list/array
        the list of orders for the same measurment as the absPos value
    @param index : list/array
        the list of indexes for the same measurment as the abPos value
    """
    return (order[absPos], index[absPos])


def order_index_selector(orders, indices, order, index, ref):
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


def group_by_numind(num, data_loc, ref):
    """
    Bins indices by the specified integer value for the argument num
    @param num : int
        number of entries per bin
    @param data_loc : str / array
        file path for data to chunk. If one wants to insert a specific array instead, put the array name in for data_loc instead.
    @param ref : int
        which measurment to use as the reference date for pulling indices. Sometimes smaller magnitude integers have nan-problems, since
        some of the earlier data on the spectrogrpahs can be a bit messy.
    @returns list_ind
        Returns a list of lists where indices are binned by the specified integer value in the num argument
    """
    # checking if the user passed a filepath into the function or a specific array
    if type(data_loc) == str:
        indices = np.array(
            hpy.File(data_loc + "All_index.hdf5", "r")["dat"][:]
        )  # points at specific indices
    else:
        indices = data_loc
    # creating a list of lists where the indices will be binned by the amount specified
    list_ind = []

    # the variable sub_list index is where the indices will be binned. list_ind is a list where all the specific bins of indices are stored
    for i in range(len(indices[ref])):
        if i % num == 0 and i != 0:
            sublist_ind = []
            for i in range(i - num, i):
                sublist_ind.append(i)
            list_ind.append(sublist_ind)

    return list_ind


def group_by_numwl(num, data_loc, ref):
    """
    Bins wavelengths by the specified integer value for the argument num
    @param num : int
        number of entries per bin
    @param data_loc : str / array
        file path for data to chunk. If one wants to insert a specific array instead, put the array name in for data_loc instead.
    @param ref : int
        which measurment to use as the reference date for pulling wavelengths. Sometimes smaller magnitude integers have nan-problems, since
        some of the earlier data on the spectrogrpahs can be a bit messy.
    @returns list_wavl
        Returns a list of lists where wavelengths are binned by the specified integer value in the num argument
    """
    # checking if the user passed a filepath into the function or a specific array
    if type(data_loc) == str:
        cenM_ar = np.array(
            hpy.File(data_loc + "All_centroidWl.hdf5", "r")["dat"][:]
        )  # points at specific location
    else:
        cenM_ar = data_loc
    # creating a list of lists where the wavlengths will be binned by the amount specified
    list_wavl = []

    # the variable sublist_wavl is where the wavelengths will be binned. list_wavl is a list where all the specific bins of wavelengths are stored
    for i in range(len(cenM_ar[ref])):
        if i % num == 0 and i != 0:
            sublist_wavl = []
            for i in range(i - num, i):
                sublist_wavl.append(cenM_ar[ref][i])
            list_wavl.append(sublist_wavl)

    return list_wavl
    

def line2Slope(line,time):
    """
    function to return the slope of a passed line
    @param line : list/array
        the intensity values for each point
    @param time : list/array
        the time values for each point
    @return a : float
        best linear fit slope for the provided data
    """
    a,b = np.polyfit(time,line,1)
    return(a)


def getAllSlopes(lines, times):
    """
    Function to fit slopes to all given lines.
    @param lines : list of 1D arrays
        Array of arrays containing intensity values for each point.
    @param times : list of 1D or 2D arrays
        Array of arrays containing time values for each point.
    @return out: list
        All of the slope values in the same order as the lines were given.
    """
    out = []
    times = np.array(times, dtype=float)  # Convert to NumPy array and ensure float type
    if times.ndim == 1:
        # If `times` is 1D, assume the same time array applies to all lines
        for i in range(len(lines)):
            out.append(line2Slope(lines[i], times))
    elif times.ndim == 2:
        # If `times` is 2D, ensure the shape matches `lines`
        for i in range(len(lines)):
            out.append(line2Slope(lines[i], times[i]))
    else:
        raise ValueError("`times` must be either 1D or 2D array.")

    return out
    

def group_by_numpix(num, data_loc, ref):
    """
    Bins wavelengths by the specified integer value for the argument num
    @param num : int
        number of entries per bin
    @param data_loc : str / array
        file path for data to chunk. If one wants to insert a specific array instead, put the array name in for data_loc instead.
    @param ref : int
        which measurment to use as the reference date for pulling wavelengths. Sometimes smaller magnitude integers have nan-problems, since
        some of the earlier data on the spectrogrpahs can be a bit messy.
    @returns list_wavl
        Returns a list of lists where wavelengths are binned by the specified integer value in the num argument
    """
    # checking if the user passed a filepath into the function or a specific array
    if type(data_loc) == str:
        cenM_pix = np.array(
            hpy.File(data_loc + "All_centroidPix.hdf5", "r")["dat"][:]
        )  # points at specific cenM_pix
    else:
        cenM_pix = data_loc
    # creating a list of lists where the pixels will be binned by the amount specified
    list_pix = []

    # the variable sublist_pix is where the pixels will be binned. list_pix is a list where all the specific bins of pixels are stored
    for i in range(len(cenM_pix[ref])):
        if i % num == 0 and i != 0:
            sublist_pix = []
            for i in range(i - num, i):
                sublist_pix.append(cenM_pix[ref][i])
            list_pix.append(sublist_pix)

    return list_pix


"""order Stuff"""

def adaptive_piecewise_linear_approximation(x, y, tolerance=1500):
    x = np.asarray(x) 
    y = np.asarray(y) 
    n = len(x)
    approx_x = [x[0]]
    approx_y = [y[0]]
    idx_start = 0  # Starting index of the current segment

    i = 1

    while i < n:
        # Fit a linear model to the current segment
        slope, intercept = np.polyfit(x[idx_start:i+1], y[idx_start:i+1], 1)

        # Calculate cumulative sum of absolute deviations from the linear fit
        cumulative_deviation = np.abs(y[idx_start:i+1] - (slope * x[idx_start:i+1] + intercept)).sum()

        # If cumulative deviation exceeds tolerance, end the segment
        if cumulative_deviation > tolerance:
            # Add the end of the segment
            approx_x.append(x[i-1])
            approx_y.append(y[i-1])
            
            # Move to the next segment
            idx_start = i - 1
        
        i += 1

    # Add the last point
    approx_x.append(x[-1])
    approx_y.append(y[-1])

    return approx_x, approx_y

def getSlopes(binnedData,Times,tolerance,doGraph=0):
    """
    calculated the slopes of each bin with respect to approximate linear eras of time. Can return slopes or plot them against bin number by era. 
    @param binnedData : list/array 
        the data from getVels that should be processed into linear eras
    @param Times : list/array
        the time array of timestamps in seconds as ints
    @param tolerance : int
        what tolerance is allowable for adaptive eras function
    @param doGraph : bol
        whether a graph should be drawn for the slopes
    """
    
    binnedData = np.concatenate(binnedData)

    import copy

    binnedBa = copy.deepcopy(binnedData)

    for i in range(len(binnedData)):
        binnedData[i] = boxcar_median(binnedBa[i],500)
        
    eras = adaptive_piecewise_linear_approximation(Times,binnedData[50],tolerance)

    fig = plt.figure()
    plt.plot(Times,binnedData[50])
    for i in eras[0]:
        plt.axvline(Times[int(i)])

    slopes = []
    for q in range(len(binnedData)):
        slopeTargets = []
        slopeTimes = []
        for i in range(1,len(eras[0])):
            slopeTargets.append(binnedData[q][int(eras[0][i-1]):int(eras[0][i])])
            slopeTimes.append(Times[int(eras[0][i-1]):int(eras[0][i])])
        soubl = []
        for i in range(len(slopeTargets)):
            soubl.append(line2Slope(slopeTargets[i],slopeTimes[i]))
        slopes.append(soubl)
        
    slopes = np.array(slopes).T[:-1]

    if(doGraph):
        fib = plt.figure()
        plt.title("slopes of drift per linear era")
        plt.xlabel("wavelength(bin number)")
        plt.ylabel("ervs drift m/s^2")
        plt.grid(1)
        for i in range(len(slopes)):
            plt.plot(boxcar_median(slopes[i],10),'.',label=i)
        plt.legend()
    return slopes


def getOrder(orders, num, bins, data):
    """
    Helper funtion to groupByOrder to chunk one order
    @param orders : array/list
        the full list of orders as pulled from the hdf5
    @param num : int
        which order to analyze
    @param bins : int
        how many bins to chop into
    @param data :
        a data file to split by the index values calculated
    @returns out : np.array
        a 2d array with an entry for each bin that contains all the data chuncked into it
    """
    mask = orders == num
    temp = []
    for i in range(len(data)):
        temp.append(data[i][mask[i]])
    filtered_data = np.array(temp)
    sizeOfBin = filtered_data.shape[1] // bins  # Size of each full bin

    trimmed_data_length = sizeOfBin * bins  # Length to trim to
    trimmed_data = filtered_data[:, :trimmed_data_length]

    out = np.split(trimmed_data, bins, axis=1)
    return np.array(out)


def groupByOrder(orders, data, bins):
    """
    Splits the set of given data into the number of bins asked for according to order.
    @param orders : list/array
        the list of orders from a hdf5 file
    @param data : list/array
        the data to chunk as pulled from hdf5 file
    @param bins : int
        how many bins you want returned
    @returns listOfBins : list
        the split data in the format [(orders),(bins),(measurments by date),(pixel/index)]
        i.e listOfBins[15,2,5,6] returns the 16th order, 3rd bin, 6th measruement, 7th index value of whatever was passed by data.
        The final dimension will likely be irregular across bins as there is no padding or cutting in place
    """
    listOfBins = []

    for q in range(
        abs(int(orders[0][0]) - int(orders[0][-1])+1)
    ):  # calculates first minus last order number to know how many orders to process
        print(
            f"Processing order {int(q+1)}/{(abs(int(orders[0][0]) - int(orders[0][-1]))+1)}",
            end="\r",
        )  # prints an estimated progress report based on the above calculation
        order_number = (
            q + orders[0][0] + 1
        )  # adjusts for order offsets like NEID Etalon starting at order 26

        dat = getOrder(orders, order_number, bins, data)
        newLis = []

        for i in range(
            len(dat)
        ):  # reformats the binned measurments and adds them to a constant list
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


def groupByOrderMeds(orders, data, bins, combine_fnc=np.nanmedian):
    """
    Splits the set of given data into the number of bins asked for according to order and takes the median within each bin per measurment.
    @param orders : list/array
        the list of orders from a hdf5 file
    @param data : list/array
        the data to chunk as pulled from hdf5 file
    @param bins : int
        how many bins you want returned
    @returns listOfBins : array
        the split data in the format [(orders),(bins),(median value of that bin per measurment)]
        i.e listOfBins[15,2,5] returns the 16th order, 3rd bin, 6th measurment's median value of whatever was passed by data.
        The final output can now be an array as shape is regular
    """
    listOfBins = []

    for q in range(abs(int(orders[0][0]) - int(orders[0][-1]))+1):
        print(
            f"Processing order {int(q+1)}/{(abs(int(orders[0][0]) - int(orders[0][-1]))+1)}",
            end="\r",
        )
        order_number = q + orders[0][0]

        dat = getOrder(orders, order_number, bins, data)
        newLis = []

        for i in range(len(dat)):
            subl = []
            for j in range(len(dat[i])):
                subl.append(
                    combine_fnc(dat[i][j])
                )  # same function as groupByOrder, but this line reduces a dimensions with medians
            newLis.append(subl)
        newList = np.array(newLis)

        listOfEras = []
        for i in range(bins):
            listOfEras.append(newList[i])
        listOfBins.append(listOfEras)

    return np.array(listOfBins)


"""Velocity Function"""


def wl2vel(wls, ref=50,ruler = None):
    """
    converts wavelength values to velocities by a reference date
    @param wls : list/array
        the wavelengths to convert
    @param ref : int (default 50)
        the reference date to convert to velocities from (getting change in wl)
    @param ruler : list/array
        optional override to the reference row used to calulate velocities, might be useful for 
        working with differences or processed wavelength data where you are looking to compare something
        other than changes since a reference date
    @returns out : list
        the list of velocities
    """
    reference_row = wls[ref].copy()
    if(ruler):
        reference_row = ruler
    out = []
    for i in range(len(wls)):
        out.append((np.subtract(wls[i], reference_row) / reference_row) * 299792458)
    return out


def getRefInds(data, num, ind):
    """
    helper function for binning a single measurment by given value
    @param data : array
        the data to binned from hdf5 file
    @param num : int
        how many pixels per bin
    @param ind :
        which measurment index to calculate
    """
    chop = data[ind][: -(len(data[ind]) % num)]
    return np.split(chop, len(chop) // num)


def getAllInds(data, num):
    """
    bins all measurments in a hdf5 file to an array
    @param data : array
        the full 2d array to bin
    @param num : int
        how many pixels per bin
    """
    out = []
    for i in range(len(data)):
        print(f"Processing {int(i+1)}/{(len(data))}", end="\r")
        out.append(getRefInds(data, num, i))
    return np.array(out)


def Fast_get_medians_in_buckets(ords, num, method=np.median):
    """
    can be manually swapped in if you are using numpy functions with axis arguments (i.e. np.median/mean/nanmedian etc.)
    works significantly faster due to unraveling, c, and dark magic.

    same docs as getRedInds
    """
    # Ensure each row's length is divisible by num
    chopped = ords[:, : (ords.shape[1] // num) * num]

    # Reshape to split each row into num-sized buckets
    reshaped = chopped.reshape(chopped.shape[0], -1, num)

    # Calculate medians along the last axis (within each bucket)
    medians = method(reshaped, axis=2)

    return medians


def get_medians_in_buckets(ords, num, agg_func=np.median):

    aggregated_values = []
    count = 0
    for row in ords:
        count += 1
        print(f"Processing {int(count+1)}/{(len(ords))}", end="\r")
        # Ensure each row's length is divisible by num
        chopped = row[: len(row) // num * num]

        # Split row into num-sized buckets
        buckets = [chopped[i : i + num] for i in range(0, len(chopped), num)]

        # Calculate aggregated value within each bucket
        bucket_aggregates = [agg_func(bucket) for bucket in buckets]

        aggregated_values.append(bucket_aggregates)

    return aggregated_values


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
    """ ""

    out = []
    data1 = np.array(data, dtype=np.float64)
    indices = np.array(indices, dtype="int")
    for i in range(len(indices)):
        out.append(data1[ref][indices[i]])
    return out


def abs_ind_vels(data, abs_ind, ref):
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
        dif = np.subtract(data[i][abs_ind], data[ref][abs_ind])
        div = np.divide(dif, data[ref][abs_ind])
        vel = np.multiply(div, c)
        velocities.append(vel)
    return velocities


def getVels(
    wavelengths, orders, bins, ordVsInd, ref=0, combine=1, combineMethod=np.nanmedian
):
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
    @param ref : int
        what index to use as reference measurment for veloicities
    @param combine : bol
        0 = return raw data | 1 = return median of each bin per measurment
    @param combineMethod : function
        the method by which to combine the lowest level of data (either order bin or index bin)
    @returns out : list
        the sliced data in terms of velocities
    """
    if not ordVsInd:
        if combine:
            return groupByOrderMeds(orders, wavelengths, bins, combineMethod)
        return groupByOrder(orders, wavelengths, bins)
    else:
        if combine:
            if (
                combineMethod.__module__ == np.__name__
            ):  # if the aggregation function is found in the numpy packadge
                print("using fast binning")
                return np.array(
                    Fast_get_medians_in_buckets(wavelengths, bins, combineMethod)
                ).T
            return np.array(get_medians_in_buckets(wavelengths, bins, combineMethod)).T
        return np.array((getAllInds(wavelengths, bins))).T

def color_gradient(vels, wavl, order, time, title, xlabel, ylabel, window_size):

    """
    Function used for creating color gradient plots of radial velocity drifts.
    @ param vels : array 
    Array of velocities that should come from arrays 
    @ param wavl : array
    Array of wavelengths that comes from  
    @ param order : array 
    array of orders for the desired dataset. Used for 
    @ param time : array 
    insert a time array for the selected data. Can be in any form except datetime values 
    @ param title : str 
    Title of the plot as indicated by the user.
    @param xlabel : str 
    Title of the x-axis as indicated by the user.
    @param ylabel : str 
    Title of the y-axis as indicated by the user. 
    @param window_size : int 
    Size of the window for taking a boxcar median of the data
    @return Gradient Plots 
    """

    #list used for storing the target velocities
    target_vels = [] 

    #for the velocities in vels, taking the boxcar median and appending it to the target_vels list
    for i in vels:
        for j in i:
            target_vels.append(boxcar_median(j,window_size))

    #turning target_vels list into an array 
    target_vels = np.array(target_vels)

    #binning the wavelengths into desired binsize
    wavl = getVels(wavl, order, 4, 0)
    median_wavl = [] 
    
    for i in range(len(wavl)):
        for j in range(4):
            med = np.nanmedian(wavl[i][j])
            median_wavl.append(med)

    # Making the meshgrid of x and y values so a gradient can be plotted
    x = time
    y = median_wavl
    X, Y = np.meshgrid(x, y)
    
    #Creating the contour plot with the color gradient 
    plt.figure(figsize=(8, 6))
    contourf_plot = plt.contourf(X, Y, target_vels, levels=20, cmap='viridis', alpha=0.7)
    
    #Plots the contour lines and puts it into the contour plot with the gradient 
    contour_lines = plt.contour(X, Y, target_vels, levels=10, cmap='viridis' , linewidths=1)
  
    plt.colorbar(contourf_plot, label='Velocities (m/s)')
    # Show plot

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

    return 

"""CSV Function"""


def out_csv(file_path, out):
    """
    Outputs velocities data into a csv file.
    @param: file_path : str
        Note: the file path cannot just be a directory it should include the desired file directory with the intended name for the file at the end of
        the directory.
    @param out : list/array
        data to store in csv

    """

    data = {"Velocities": out}

    df = pd.DataFrame(data)
    df.to_csv(file_path, index=False)

    return
