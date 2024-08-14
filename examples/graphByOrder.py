#import the project library
import lib

#import everything else
import h5py as hpy
import matplotlib.pyplot as plt

#open data files that have been pre-processed
wavl = hpy.File("./All_centroidWl.hdf5",'r')['dat'][:]
ords = hpy.File("./All_order.hdf5",'r')['dat'][:]

#(optional) convert wavelegnths to equivelent velocities from a reference date
wavl = lib.wl2vel(wavl)

#process data into per-order wavelength bins
ordVel = lib.getVels(wavl,ords,5,0)
arr = np.array(ordVel)

#(optional) smooth data with a moving median function
def boxcar_median(data, window_size):
    half_window = window_size // 2
    smoothed_data = np.zeros(len(data), dtype=np.float64)
    
    for i in range(len(data)):
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        smoothed_data[i] = np.median(data[start:end])
    
    return smoothed_data


#plot the data
fig = plt.figure()
for i in range(5):
    #boxcar_median may be removed here
    plt.plot(boxcar_median(arr[50][i],100),label=i+1)
plt.title("51st order binned to 5 using github library")
plt.legend()
plt.xlabel("measurment number")
plt.ylabel("ervs")
