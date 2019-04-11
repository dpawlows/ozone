'''Plot output vs time'''
import glob
import pdb
from srcOutput import readOutput
from matplotlib import pyplot as pp
import numpy as np
import datetime
import matplotlib.dates as mdates
import os

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

f = open(files[0],'r')

altitude = []
i=0
started = False
for line in f:
    if started:
        t = line.split()
        altitude.append(float(t[0]))
        print('{}: {}'.format(i,altitude[i]))
        i+=1

    if line[0:4] == "#Alt":
        vars = line.strip("#").strip("\n").split("\t")
        started = True

iAlt = int(input("Enter which alt to plot: "))
f.close()
for i in range(1,len(vars)):
    #skip altitude as an option
    print("{}: {}".format(i,vars[i]))
iVar = int(input("Enter which variable to plot: (-1 for all)"))

time = []
data = []
dPlanet = []
for file in files:
    temp,header = readOutput.readOutputFile(file)
    dPlanet.append(float(header[-1].split("=")[1]))
    if iVar < 0:
        data.append(temp[iAlt,:])
    else:
        data.append(temp[iAlt,iVar])
    pos = file.find(".dat") - 12
    days = int(file[pos:pos+5])
    hours = int(file[pos+6:pos+8])
    mins = int(file[pos+8:pos+10])
    secs = int(file[pos+10:pos+12])


    time.append(datetime.timedelta(days=days,hours=hours,\
    minutes=mins,seconds=secs))

data = np.array(data)

ymin = float(input("Enter minimum to plot (0 for auto):"))
ymax = float(input("Enter maximum to plot (0 for auto):"))

time = [t.days*24+t.seconds/3600. for t in time]
plotDistance = True
fig=pp.figure()
ax = fig.add_subplot(221)
if iVar < 0:
    for i in range(1,len(vars)):
        pp.plot(time,data[:,i],label=vars[i])
    pp.ylabel(vars[1][0])

else:
    pp.plot(time,data)
    pp.ylabel(vars[iVar])
if ymin != 0 and ymax != 0:
    pp.ylim([ymin,ymax])
# pp.ylim([2.17e12,2.18e12])
if plotDistance:
    ax = fig.add_subplot(223)
    pp.plot(time,dPlanet)
    pp.ylabel("Planet Distance (AU)")
    pp.xlabel("Time (hours)")
else:
    pp.tight_layout()
    pp.xlabel("Time (hours)")

pp.savefig(datadir+'/plot.png')
