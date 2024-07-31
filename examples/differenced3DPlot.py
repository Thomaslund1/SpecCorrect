import h5py as hpy
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import lib
%matplotlib widget



sciWaves = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5",'r')['dat'][:]
ords = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5",'r')['dat'][:]
sciTimes = np.array(hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_time.hdf5",'r')['dat'][:],dtype=str)
calWaves = hpy.File("/data1/linefits/outputsEtalonNEID/SciCurrent/All_centroidWl.hdf5",'r')['dat'][:]
calTimes = np.array(hpy.File("/data1/linefits/outputsEtalonNEID/SciCurrent/All_time.hdf5",'r')['dat'][:],dtype=str)

sciWaves = lib.interpolate(sciWaves)
calWaves = lib.interpolate(calWaves)

sVEL = []
for i in range(len(sciWaves)-2):
    sVEL.append((np.subtract(sciWaves[-1],sciWaves[i])/sciWaves[-1])*(3*10**8))
cVEL = []
for i in range(len(calWaves)-2):
    cVEL.append((np.subtract(calWaves[-1],calWaves[i])/calWaves[-1])*(3*10**8))


wlMeds = []
for i in sVEL:
    wlMeds.append(np.nanmedian(i))

calWlMeds = []
for i in cVEL:
    calWlMeds.append(np.nanmedian(i))

times = []
for i in range(len(sciTimes)):
    times.append(datetime.fromisoformat(str(sciTimes[i][:8] + str("T") + sciTimes[i][8:])).timestamp())

def boxcar_median(data, window_size):
    half_window = window_size // 2
    smoothed_data = np.zeros(len(data), dtype=np.float64)
    
    for i in range(len(data)):
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        smoothed_data[i] = np.median(data[start:end])
    
    return smoothed_data

diffs = np.subtract(sVEL,cVEL[:])


diffs = lib.getVels(diffs,ords,5,0)

diffs2 = []
for i in diffs:
    for j in i:
        diffs2.append(lib.boxcar_median(j,700))
diffs2 = np.array(diffs2)
print(diffs2.shape)



diffs2 = np.array(diffs2,dtype='float')

diffs2 = lib.interpolate(diffs2)


y = np.arange(diffs2.shape[0])  # This will be [0, 1, 2, ..., 399]
x = times[:-2]# Time array

# Create meshgrid for x (time) and y (rows of the data)
x, y = np.meshgrid(x, y)

# Create the figure and 3D axes
fig = plt.figure(figsize=(14, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
ax.plot_surface(x, y, diffs2, cmap='viridis', edgecolor='none', rstride=8, cstride=8)

# Set labels
ax.set_xlabel('Time (s))')
ax.set_ylabel('wavelength (bin#)')
ax.set_zlabel('ERVS')
ax.set_title('early vs late reference file differences')

# Show the plot
plt.show()

