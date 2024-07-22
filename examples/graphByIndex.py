import lib
import matplotlib.pyplot as plt
import h5py as hpy

#grab data
wavl = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5",'r')['dat'][:]
ords = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5",'r')['dat'][:]

#process data by index
ordVel = getVels(wavl,ords,15,1)

#graph a few bins
fig = plt.figure()
for i in range(5):
    plt.plot(boxcar_median(arr[1500+(i*10)],100),label=i)
plt.xlabel("measurment number")
plt.ylabel("ervs")
