'''Plot output vs time'''
import glob
import pdb
from srcOutput import readOutput
from matplotlib import pyplot as pp
import datetime
import matplotlib.dates as mdates
import os

newfiles = glob.glob("data/*dat")
files = sorted( newfiles, key = lambda file: os.path.getctime(file))
outputPlotDir = ("data/plots/")
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

fig=pp.figure()
ax = fig.add_subplot(221)
pp.tight_layout()

for iSpecies in range(1,len(vars)):
    ax.semilogx(data[:,iSpecies],data[:,0],
        label=vars[iSpecies])
pp.xlabel('Density (#/cm3)')
pp.ylabel('Altitude (km)')
pp.xlim([1e-9,1e22])
pp.legend(loc='upper left',frameon=False)
pp.text(.02, 1.02, '{}d'.format(time.days),
 transform=ax.transAxes)
pp.savefig('profile.png')
