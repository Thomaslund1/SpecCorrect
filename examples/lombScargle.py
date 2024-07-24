#import the library
import lib

#import everything else
import numpy as np
import h5py as hpy
import matplotlib.pyplot as plt
from datetime import datetime

#pull needed data
waves = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5",'r')['dat'][:]
ords = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5",'r')['dat'][:]
timesF = np.array(hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_time.hdf5",'r')['dat'][:],dtype='str')

#convert datetimes from iso_T strings to unix timestamps
times = []
for i in range(len(timesF)):
    times.append((datetime.fromisoformat(str(timesF[i][:8] + str("T") + timesF[i][8:]))).timestamp())

#get binned data by order 
byOrder = lib.getVels(waves,ords,4,0)

#plot the reference data
fig = plt.figure()
for i in range(4):
    plt.plot(times,byOrder[0][i],label=i)
    plt.grid(True)
    plt.legend()

#make a lomb scargle periodogram of one of ths bins with a logarithmic periods axis from 10**-3 to 10**3
lib.lSPerodogram(np.logspace(-3,3,1000),times,byOrder[0][2],1)
