"""Compare two output files to see if there are significant
differences between the two files"""
import glob
import pdb
import numpy as np
def load_src(name, fpath):
    import os, imp
    p = fpath if os.path.isabs(fpath) \
        else os.path.join(os.path.dirname(__file__), fpath)
    return imp.load_source(name, p)

load_src("readOutput","../srcOutput/readOutput.py")
from readOutput import readOutputFile

testfile = "outputtest00160_000000.dat"
files = glob.glob("output?????_??????.dat")
for i in range(len(files)):
    print("{}: {}".format(i,files[i]))
ifile = int(input("Enter which file to compare:"))

print("Comparing {} with {}".format(files[ifile],testfile))

data,vars = readOutputFile(files[ifile])
tdata,tvars = readOutputFile(testfile)

for iSpecies in range(1,len(vars)):
    itvar = tvars.index(vars[iSpecies])
    tvdata = tdata[:,itvar]
    vdata = data[:,iSpecies]
    diff = abs((tvdata-vdata))/tvdata
    maxdiff = max(diff)
    if maxdiff > .1:

        print("Test Failed!!!")
        print("Max difference of {} found at index {}, species {}"\
        .format(maxdiff,np.where(diff==maxdiff)[0],vars[iSpecies]))
        exit(1)

print("Test Passed!!!")
