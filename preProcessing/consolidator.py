import numpy as np
import glob
from datetime import datetime
import multiprocessing
import h5py as hpy
import sys


def getFiles(Path):
    return(np.array(np.sort(glob.glob(Path + "/*.npz"))))


def loadData(file_index, header_num,Files):
    """
    pulls all data from the given header and returns it as an array
    """
    dat_file = Files[file_index]
    with np.load(dat_file, allow_pickle=1) as dat:
        return dat['arr'][header_num]
        
def getDate(fileIndex,Files):
    """
    returns a list of seconds since Jan 1 2020 at 00:00
    for each given file
    """
    refDate = datetime(2020,1,1)
    dateTime = re.findall('\d{8}T\d{6}',str(file))[0]
    dateForm = "%Y%m%dT%H%M%S"
    DateObj = datetime.strptime(dateTime,dateForm)
    return(DateObj)


def process_file(file_index,header,Files):
    """
    helper function to grab wavelength and timestamp for one file
    """
    out_wl = loadData(file_index, header,Files)
    return out_wl

def timeTest(column_index,allFiles,cores,path):
    """
    Multithreaded method to process files in parallel and store results in numpy arrays.
    """
    length = len(allFiles)
    batch_size = 100
    num_batches = length // batch_size
    
    pool = multiprocessing.Pool(cores)
    
    # Use lists to collect results
    out_wls = []

    print(range(num_batches))
    
    for i in range(num_batches):
        print(f"Processing batch {i+1}/{num_batches}", end="\r")
        start_index = i * batch_size
        end_index = min((i + 1) * batch_size, length)
        batch_indices = range(start_index, end_index)
        
        # Use pool.starmap to parallelize process_file
        batch_results = pool.starmap(process_file, [(idx, column_index,allFiles) for idx in batch_indices])
        
        # Collect batch_results in lists
        for wl in batch_results:
            out_wls.append(wl)
    
    pool.close()
    pool.join()
    
    # Convert lists to numpy arrays
    out_wls = np.array(out_wls)
    
    return out_wls



def ConsolidateCol(path,allFiles,cores):
    keys = ["order","index","centroidPix","centroidWl","fwhmPix","fwhmWls","snrPeak"]
    for i in range(6
        print("Currently Crunching " + keys[i] + ". Expect the first collum to take the longest")
        savePath = path + "/All_" + keys[i] + ".hdf5"
        outFile = hpy.File(savePath,'w')
        outArr = timeTest(i,allFiles,cores,savePath)
        outFile.create_dataset("dat",data=outArr)
        outFile.close
        print("column complete!")
    timeList = []
    print("generating reference time list: ")
    for i in range(len(allFiles)):
        timeList.append(getDate(i,allFiles))
    timeList = np.array(timeList)
    savePath = path + "/All_RefTimes.hdf5"
    outFile = hpy.File(savePath,'w')
    outFile.create_dataset("dat",data=timeList)
    outFile.close()



def main():
    args = [None,None,1]
    inputs = sys.argv[1:]
    if(len(inputs) < 2):
        print("Usage: consolodate.py {inputs directory} {outputs directory} {number of cores (only helps with sufficent free memory available)}")
        return
    for i in range(len(inputs)):
        args[i] = inputs[i]
    print(args)
    for i in args:
        if i == None:
            print("Usage: consolodate.py {inputs directory} {outputs directory} {number of cores (only helps with sufficent free memory available)}")
            return
    print("Pulling Data from files")
    FilesList = getFiles(args[0])
    ConsolidateCol(args[1],FilesList,int(args[2]))

if(__name__ == '__main__'):
    main()
