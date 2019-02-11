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

    if line[0] == "#":
        vars = line.strip("#").strip("\n").split("\t")
        started = True

iAlt = int(input("Enter which alt to plot: "))
f.close()
for i in range(1,len(vars)):
    #skip altitude as an option
    print("{}: {}".format(i,vars[i]))

iVar = int(input("Enter which variable to plot: "))

time = []
data = []

for file in files:
    temp,header = readOutput.readOutputFile(file)
    data.append(temp[iAlt,iVar])
    pos = file.find(".dat") - 12
    days = int(file[pos:pos+5])
    hours = int(file[pos+6:pos+8])
    mins = int(file[pos+8:pos+10])
    secs = int(file[pos+10:pos+12])


    time.append(datetime.timedelta(days=days,hours=hours,\
    minutes=mins,seconds=secs))

time = [t.days*24+t.seconds/3600. for t in time]
fig=pp.figure()
ax = fig.add_subplot(221)
pp.tight_layout()
pp.plot(time,data)
pp.xlabel("Time (hours)")
pp.ylabel(vars[iVar])
pp.savefig('plot.png')
