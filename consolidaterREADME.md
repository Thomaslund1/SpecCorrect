Consolodater / ConsolodaterMemEdition 
=======
Jun 20 2024

Overview
-----

These programs are designed to process through the many .npz output files given by large linefits runs and packadge them into much friendlier database formats for further procecssing. This step reduces the acess time and memory usage for pulling data across many measurments by several orders of magnitude. 

Usage
----
>python Consolodater {Str: Input Directory} {Str: Output Directory} {number of cores = 1}

>python ConsolodaterMemEdition.py {Str: Input Directory} {Str: Output Directory}

The consolodater will output to 7 .hdf5 files each with one master array called 'dat' all the data will be indexed in chronological order under that header  

There is one additional output beyond the 6 columns found in the .npz files for a standard date scheme. This outputs the time of each file as an int for easy storage. The date format is YYYYMMDDHHmmSS. The same as the date string format just without the 'T' to seperate the time string. 

the Consolodater standard edition can use multithreading well, but can take very large chunks of ram so use with caution. Consolodater memory edition will chunk the files by hundreds, append to the output, then when all files are processed, combine each chunk of 100 by column.

to open an output file:
----
>import h5py as hpy

>file = hpy.File("{path}","{r/w/a/x/etc}")

>data = file['dat'][:]



to write to an output file with more date
----
>import h5py as hpy

>file = hpy.File("{path}","a")

>data = file.create_dataset("dat{x}",data={data goes here})
