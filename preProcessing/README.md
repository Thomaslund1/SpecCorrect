Consolodater / ConsolodaterMemEdition 
=======
Jun 20, 2024

Overview
-----

These programs are designed to process through the many .npz output files given by large linefits runs and package them into much friendlier database formats for further processing. This step reduces the access time and memory usage for pulling data across many measurements by several orders of magnitude. 

Usage
----
>python Consolodater {Str: Input Directory} {Str: Output Directory} {number of cores = 1}

>python ConsolodaterMemEdition.py {Str: Input Directory} {Str: Output Directory}

The consolidator will output to 7 .hdf5 files each with one master array called 'dat' all the data will be indexed in chronological order under that header  

There is one additional output beyond the 6 columns found in the .npz files for a standard date scheme. This outputs the time of each file as an int for easy storage. The date format is YYYYMMDDHHmmSS. The same as the date string format just without the 'T' to separate the time string. 

the Consolidator standard edition can use multithreading well, but can take very large chunks of ram so use with caution. Consolidator memory edition will chunk the files by hundreds, append them to the output, and then when all files are processed, combine each chunk of 100 by column.

to open an output file:
----
>import h5py as hpy

>data = hpy.File("{path}","{r/w/a/x/etc}")['dat'][:]


runTest
------
A quick test program to aid linefits runs. Checks the input and output directory for matching date strings in file names and can generate a list of discrepencies.
useful for things like checking the progress of a linefits run, reparing a failed or interrupted run, or setting target lists
