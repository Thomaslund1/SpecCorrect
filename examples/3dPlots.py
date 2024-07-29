import h5py as hpy
import numpy as np
import matplotlib.pyplot as plt
import ryanFunctions as lib

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

#convert wavelengths to velocities
ervs1 = lib.wl2vel(waves1)
ervs2 = lib.wl2vel(waves2)

ervs1 = lib.getVels(ervs1,ords1,5,0)
ervs2 = lib.getVels(ervs2,ords2,5,0)

from datetime import datetime

Times1 = []
for i in range(len(times)):
    Times1.append((datetime.fromisoformat(str(times[i][:8] + str("T") + times[i][8:]))).timestamp())
Times2 = []
for i in range(len(times1)):
    Times2.append((datetime.fromisoformat(str(times1[i][:8] + str("T") + times1[i][8:]))).timestamp())

diffs = np.subtract(ervs1,ervs2)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import h5py as hpy
%matplotlib widget


newErvs1 = []
for i in diffs:
    for j in i:
        newErvs1.append(lib.boxcar_median(j,50))
newErvs1 = np.array(newErvs1)

newTimes = []
for i in range(400):
    newTimes.append(Times1)
newTimes = np.array(newTimes)
newTimes.shape


y = np.arange(newErvs1.shape[0])  # This will be [0, 1, 2, ..., 399]
x = Times1  # Time array

# Create meshgrid for x (time) and y (rows of the data)
x, y = np.meshgrid(x, y)

# Create the figure and 3D axes
fig = plt.figure(figsize=(14, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
ax.plot_surface(x, y, newErvs1, cmap='viridis', edgecolor='none', rstride=8, cstride=8)

# Set labels
ax.set_xlabel('Time')
ax.set_ylabel('Row Index')
ax.set_zlabel('Value')
ax.set_title('3D Plot of master file differences')

# Show the plot
plt.show()
