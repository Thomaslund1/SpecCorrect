import numpy as np
import glob
import re
import sys

def getList(filePathA,filePathB):
    """
    takes two file paths and returns a chronologically sorted np array of the inputs and outputs
    """
    return(np.sort(glob.glob(filePathA + "/*.fits")),np.sort(glob.glob(filePathB + "/*.npz")))


def getDates(listOfFiles):
    """
    takes a file list and returns a list of datetime strings pulled from those file names 
    """
    dates = []
    for i in listOfFiles:
        dates.append(re.findall('\d{8}T\d{6}',str(i))[0])
    return np.sort(dates)

def checkFiles(DatesB,DatesA):
    """
    takes two lists of dates and checks each value from one list against the other. Shared items are removed. The remaining list and a list of files that didnt match 
    are returned
    """
    worked = []
    didntWork = []
    DatesA = list(DatesA)
    DatesB = list(DatesB)
    for i in range(len(DatesA)):
        if(DatesA[i] in DatesB):
            worked.append(i)
            DatesB.remove(DatesA[i])
        else:
            didntWork.append(i)
    print("Files in output folder: ",len(DatesA))
    print("Files in reference but not output folder: ",len(DatesB))
    print("Correct files in output folder: ",len(worked))
    print("Incorrect files in output folder: ",len(didntWork))
    return DatesB,didntWork

def getMissed(remainderList,inptPath):
    """
    Writes the missed files to a file according to user input
    """
    name = input("Write name of output file (i.e TargetFiles.txt )\n")
    print("collecting missing files. This is the longest part of the process hold on...")
    missed = []
    for i in remainderList:
        missed.append(glob.glob(inptPath + "/*" + i + "*.fits" ))
    with open(name,'w') as file:
        for i in missed:
            file.write(str(i)[2:-2] + '\n')
    print("targets written to file")

def main():
    args = [None,None,0]
    for i in range(len(sys.argv[1:])):
        args[i] = sys.argv[i+1]
    for i in args:
        if i == None:
            print("Usage: runTest.py {Input directory} {output directory} {generate list of targets = 0}")
            return
    args[2] = int(sys.argv[3])
    print(args)
    fileA,fileB = getList(args[0],args[1])
    A = getDates(fileA)
    B = getDates(fileB)
    remDats,noWork = checkFiles(A,B)
    if(args[2]):
        getMissed(remDats,args[0])

if __name__ == "__main__":
    main()
