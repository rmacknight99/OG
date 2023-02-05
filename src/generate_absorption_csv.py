from os.path import isfile, isdir, join, dirname, realpath, basename, exists
from os import listdir, mkdir
from argparse import ArgumentParser
import pandas as pd
import numpy as np

parser = ArgumentParser()
parser.add_argument('--data',type=str,required=True,help='file or directory containing tddft data in orca 5 format')
parser.add_argument('--ext', type=str,default=None,help='common extension used on all outfiles')
parser.add_argument('--roots',type=int,required=True,help='Number of roots or transitions used in the tddft calculations')
parser.add_argument('--outpath',type=str,default=None,help='outpath for generated csv(s)')
args=parser.parse_args()

data_path=args.data
out_path=args.outpath
extension=args.ext
roots=args.roots


# Converts a list of strings into a rootsx3 array 
def list_2_array(lst):
    array=np.empty(shape=(len(lst),3))
    for i,line in enumerate(lst):
        data=line.split()
        array[i,:]=data[:3]
    return array


# If no extension provided defaults to an out file
if extension is None:
    extension='.out'

is_directory=False
if isdir(data_path):
    files=[f for f in listdir(data_path) if isfile(join(data_path,f)) and extension in f]
    is_directory=True
elif isfile(data_path):
    files=data_path
else:
    raise Exception('Provided path for data is neither a directory or file.')

if out_path is None:
    if not is_directory:
        out_path=dirname(realpath(data_path)) 
    else:
        out_path=data_path

# Make the output directory if it does not exist
if is_directory and not out_path is None:
    if not exists(out_path):
        mkdir(out_path)


# Generate csvs 
if not is_directory:
    bname=basename(files).replace(extension,'.csv')
    out_path=join(out_path,bname)
    with open(files) as f:
        lst=f.readlines()
    for n, line in enumerate(lst):
        #IMPORTANT
        #This assumes the format of a tddft outfile from orca 5 
        #Uses a different software or version likely makes this script useless 
        if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
            spectra=lst[n+5:n+5+roots]
            spectra=list_2_array(spectra)
            break
    df=pd.DataFrame(data=spectra,columns=['State','Energy (cm-1)','Wavelength (nm)'])
    df.to_csv(path_or_buf=out_path)

else:
    for out_file in files:
        bname=basename(out_file).replace(extension,'.csv')
        out_csv=join(out_path,bname)
        with open(join(data_path,out_file)) as f:
            lst=f.readlines()
        for n, line in enumerate(lst):
            #IMPORTANT
            #This assumes the format of a tddft outfile from orca 5 
            #Uses a different software or version likely makes this script useless 
            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                spectra=lst[n+5:n+5+roots]
                spectra=list_2_array(spectra)
                break
        df=pd.DataFrame(data=spectra,columns=['State','Energy (cm-1)','Wavelength (nm)'])
        df.to_csv(path_or_buf=out_csv)






    






#          ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
