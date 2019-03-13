'''Plot output vs time'''
import glob
import pdb
from srcOutput import readOutput
from matplotlib import pyplot as pp
import datetime
import matplotlib.dates as mdates
import numpy as np
# from matplotlib import rc
import os

pp.rc('text', usetex=True)
datadir = input("Enter data directory: ")
newfiles = glob.glob(datadir+"/*dat")
filetypes = set([])
for file in newfiles:
    pos = file.rfind("/")+1
    filetypes.add(file[pos:pos+5])
filetypes = list(filetypes)
filetypes.sort()

for i in range(len(filetypes)):
    print(i,filetypes[i])
if len(filetypes) > 1:
    itype = int(input("Enter which type of file to plot: "))
    newfiles = glob.glob(datadir+"/"+filetypes[itype]+"*.dat")

files = sorted( newfiles, key = lambda file: os.path.getctime(file))
nfiles = len(files)

if nfiles < 1:
    print("No files found!")
    exit(1)

for i in range(nfiles):
    print('{}: {}'.format(i,files[i]))
ifile = int(input("Enter which file to plot: "))

data,vars = readOutput.readOutputFile(files[ifile])

pos = files[ifile].find(".dat") - 12
days = int(files[ifile][pos:pos+5])
hours = int(files[ifile][pos+6:pos+8])
mins = int(files[ifile][pos+8:pos+10])
secs = int(files[ifile][pos+10:pos+12])
time = datetime.timedelta(days=days,hours=hours,\
minutes=mins,seconds=secs)

for i in range(len(vars)):
    print(i,vars[i])
nvarstoplot = int(input("Enter how many vars to plot: "))
pvars = []
for i in range(nvarstoplot):
    pvars.append(int(input("Enter var "+str(i)+": ")))

fig=pp.figure()
ax = fig.add_subplot(221)
pp.tight_layout()
# tempvar = [1 for i in vars if i[0] == '[']

for iSpecies in pvars:
    ax.semilogx(data[:,iSpecies],data[:,0],
        label=vars[iSpecies])
pp.xlabel('Sources and losses (cm$^{-3}$s$^{-1}$)')
pp.ylabel('Altitude (km)')
xmin = float(input("Enter minimum to plot (0 for auto):"))
xmax = float(input("Enter maximum to plot (0 for auto):"))

if xmin != 0 and xmax != 0:
    pp.xlim([xmin,xmax])
pp.legend(loc='upper left',frameon=False)
pp.text(.02, 1.02, '{}d'.format(time.days),
 transform=ax.transAxes)
pp.savefig(datadir+'/profile.png')
