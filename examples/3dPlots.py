import h5py as hpy
import numpy as np
import matplotlib.pyplot as plt
import ryanFunctions as lib
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D


#collect all needed files
waves1 = hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5','r')['dat'][:]
waves2 = hpy.File('/data1/linefits/outputsEtalonNEID/newMasterConfig/consolodated_current/All_centroidWl.hdf5','r')['dat'][:]
ords1  = hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5','r')['dat'][:]
ords2  = hpy.File('/data1/linefits/outputsEtalonNEID/newMasterConfig/consolodated_current/All_order.hdf5','r')['dat'][:]
times = np.array(hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_time.hdf5','r')['dat'][:],dtype=str)
times1 = np.array(hpy.File('/data1/linefits/outputsEtalonNEID/newMasterConfig/consolodated_current/All_time.hdf5','r')['dat'][:],dtype=str)

#interpolate to remove nans(optional)
waves1 = lib.interpolate(waves1)
waves2 = lib.interpolate(waves2)
ords1 = lib.interpolate(ords1)
ords2 = lib.interpolate(ords2)

#convert wavelengths to velocities (comparing relative drift differences between fibers instead of absolute drifts)
#to see absolute drift differences one might subtract wavelengths, then run that array through wl2vel with a ruler array of wavelengths from one of the initial wavelength arrays
ervs1 = lib.wl2vel(waves1)
ervs2 = lib.wl2vel(waves2)

#bin the velocities
ervs1 = lib.getVels(ervs1,ords1,5,0)
ervs2 = lib.getVels(ervs2,ords2,5,0)

#convert times to unix timestamps
Times1 = []
for i in range(len(times)):
    Times1.append((datetime.fromisoformat(str(times[i][:8] + str("T") + times[i][8:]))).timestamp())

#subtract the binned velocities to  see how they differ
diffs = np.subtract(ervs1,ervs2)

#flatten the array and take a boxcar median for outliers
newErvs1 = []
for i in diffs:
    for j in i:
        newErvs1.append(lib.boxcar_median(j,50))
newErvs1 = np.array(newErvs1)

#duplicate the times array to match the dimensions of our data array
newTimes = []
for i in range(400):
    newTimes.append(Times1)
newTimes = np.array(newTimes)
newTimes.shape

#fancy 3d plots
y = np.arange(newErvs1.shape[0])  # create one axes for the length of our data array (these will be used for order bins in this case)
x = Times1 # next axis is for time
x, y = np.meshgrid(x, y) # turn to meshgrid for plotting
fig = plt.figure(figsize=(14, 7)) # make a figure
ax = fig.add_subplot(111, projection='3d') # matplotlib 3d schenangians 
ax.plot_surface(x, y, newErvs1, cmap='viridis', edgecolor='none', rstride=8, cstride=8) # plot the x and y axes from earlier against our data array

#formatting
ax.set_xlabel('Time')
ax.set_ylabel('Row Index')
ax.set_zlabel('Value')
ax.set_title('3D Plot of master file differences')
plt.show()
