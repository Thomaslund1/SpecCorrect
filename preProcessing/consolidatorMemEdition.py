import numpy as np
import glob
from datetime import datetime
import multiprocessing
import h5py as hpy
import sys
import re

def process_chunk(chunk_files,keys,outDir,iteration):
    """
    Helper function to process one chunk of files passed by process_files_in_chunks.
    """
    num_files = len(chunk_files)
    out_wls = {key: [] for key in keys}  # Dictionary to store results for each key
    
    for file_path in chunk_files:
        # Process each file and collect results for all keys
        # Example: Assuming process_file returns a dictionary with keys as specified
        file_data = loadData(file_path)
        for key in keys:
            out_wls[key].append(file_data[key])
    
    # Convert lists to numpy arrays
    out_arrs = {key: np.array(value) for key, value in out_wls.items()}
    
    # Update HDF5 files for each key
    for key in keys:
        save_path = str(outDir + "/Temp_" + key + ".hdf5") 
        print(save_path)
        try:
            doesExist = hpy.File(save_path, 'a')
            doesExist.close()
        except:
            doesExist = hpy.File(save_path,'w')
            doesExist.close()
        with hpy.File(save_path, 'a') as outFile:
            outFile.create_dataset("dat" + iteration, data=out_arrs[key])

def process_files_in_chunks(allFiles,outDir,chunk_size=100):
    """
    Process all files in chunks and update HDF5 files for each column.
    """
    keys = ["order", "index", "centroidPix", "centroidWl", "fwhmPix", "fwhmWls", "snrPeak","time"]
    length = len(allFiles)
    num_chunks = (length + chunk_size - 1) // chunk_size  # Calculate number of chunks
    
    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, length)
        chunk_files = allFiles[start_idx:end_idx]
        
        print(f"Processing chunk {chunk_idx + 1}/{num_chunks}")
        process_chunk(chunk_files, keys,outDir,str(chunk_idx))
        print(f"Chunk {chunk_idx + 1} complete")


# Example function to simulate processing a file
def loadData(file_path):
    """
    pulls all data from the given header and returns it as an array
    """
    dat_file = file_path
    with np.load(dat_file, allow_pickle=1) as dat:
        time = getDate(file_path)
            
        return {"order": dat['arr'][0],
            "index": dat['arr'][1],
            "centroidPix": dat['arr'][2],
            "centroidWl": dat['arr'][3],
            "fwhmPix": dat['arr'][4],
            "fwhmWls": dat['arr'][5],
            "snrPeak": dat['arr'][6],
            "time": time}

def getFiles(Path):
    return(np.array(np.sort(glob.glob(Path + "/*.npz"))))

def getDate(file):
    """
    returns a list of seconds since Jan 1 2020 at 00:00
    for each given file
    """
    dateTime = re.findall('\d{8}T\d{6}',str(file))[0]
    newDateTime = dateTime[:8] + dateTime[9:]
    return(int(newDateTime))

def PS(outDir):
    """
    Post Sript to turn chunks of data into single columns. Optional but useful for data where one optimized collum at a time can fit in memory.
    If disabled all data will be stored in the Temp_{header}.hdf5 files by group of 100 measurment
    """
    print("Cleaning up some data, this may take a while...")
    keys = ["order", "index", "centroidPix", "centroidWl", "fwhmPix", "fwhmWls", "snrPeak", "time"]
    for key in keys:
        save_path = outDir + "/Temp_" + key + ".hdf5"
        correction_path = outDir + "/All_" + key + ".hdf5"
        with hpy.File(save_path, 'r+') as file:
            combined_data = []
            for i in range(len(file.keys())):  # Iterate over all keys
                dataset_name = str('dat' + str(i))
                if dataset_name in file:
                    combined_data.append(file[dataset_name][:])
                
            # Concatenate all datasets into a single numpy array
            if combined_data:
                ar = np.concatenate(combined_data)

                # Create or overwrite the "dat" dataset with the combined data
                with hpy.File(save_path, 'w') as finalOut:
                    finalOut.create_dataset("dat", data=ar)
                    

    print("Done")
    return
    


# Example usage:
def main():
    args = [None,None,1]
    inputs = sys.argv[1:]
    if(len(inputs) < 1):
        print("Usage: consolodate.py {inputs directory} {outputs directory}")
        return
    for i in range(len(inputs)):
        args[i] = inputs[i]
    print(args)
    for i in args:
        if i == None:
            print("Usage: consolodate.py {inputs directory} {outputs directory}")
            return
    '''
    MakeDir(args[1])
    '''
    allFiles = getFiles(args[0])  # Replace with your list of files
    process_files_in_chunks(allFiles,args[1], chunk_size=100)

    #may be commented out
    PS(args[1])

if (__name__ == "__main__"):
    main()
