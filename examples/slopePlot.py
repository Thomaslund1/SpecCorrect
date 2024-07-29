import lib
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import h5py as hpy

#pull data
waves = hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5','r')['dat'][:]
ords = hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5','r')['dat'][:]
times = timesCal = np.array(hpy.File('/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_time.hdf5','r')['dat'][:],dtype=str)

#process and bin
waves = lib.interpolate(waves)
waves = lib.getVels(waves,ords,5,0)

#convert datetimes to timestamps
Times1 = []
for i in range(len(times)):
    Times1.append((datetime.fromisoformat(str(times[i][:8] + str("T") + times[i][8:]))).timestamp())

#make everything an array
Times1 = np.array(Times1)
waves = np.array(waves)

#reduce a dimension of the wavelength data
newWaves = []
for i in waves:
    for j in i:
        newWaves.append(j)
newWaves = np.array(Times1)

#calculate all slopes for the orders over time
slopes = lib.getAllSlopes(newWaves,newTimes)

#plot slopes against bin number
fig = plt.figure()
plt.plot(slopes)
