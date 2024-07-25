#standard operations
import lib
waves = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_centroidWl.hdf5",'r')['dat'][:]
ords = hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_order.hdf5",'r')['dat'][:]
timesF = np.array(hpy.File("/data1/linefits/outputsEtalonNEID/consolidatedData_current/All_time.hdf5",'r')['dat'][:],dtype='str')
times = []
for i in range(len(timesF)):
    times.append((datetime.fromisoformat(str(timesF[i][:8] + str("T") + timesF[i][8:]))).timestamp())

#binning by wave
byWave = getVels(waves,ords,50,1)

#as the structure is [bins][measurments] if we flip the array sideways we get the full readout from one measurment for each bin
#using the first measure as a spectral ruler, we can compare the first to last measurment change in wavelength over wavelength
ruler = byWave.T[0]
fig = plt.figure()
plt.plot(ruler,np.subtract(byWave.T[0],byWave.T[-1]))


#same example for order bins

byOrder = lib.getVels(waves,ords,5,0)

ruler = byOrder_np[:, :, 0].flatten()
