from itertools import islice
import numpy as np
from matplotlib import pyplot as pp
import textwrap as tw

#Variables I am tracking
pressure = []
starTemp = []
atmTemp = []
maxOzone = []
altMax = []
totalColumn = []

#Pull pressures into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
     for line in islice(f, 1, None, 6):
        pressure.append(float(line))

#Pull star temperatures into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
    for line in islice(f, 2, None, 6):
        starTemp.append(int(line))

#Pull average atmosphere temperatuers into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
    for line in islice(f, 3, None, 6):
        atmTemp.append(float(line))

#Pull maximum ozone into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
    for line in islice(f, 4, None, 6):
        maxOzone.append(float(line))

#Pull altitude of maximum ozone into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
     for line in islice(f, 5, None, 6):
         altMax.append(float(line))

#Pull total columns into array
with open('runMultipleEOMdata/runMultipleLog.txt') as f:
     for line in islice(f, 6, None, 6):
         totalColumn.append(float(line))

#Plot data
fig, ax1 = pp.subplots(1, sharex=True)
scatter1 = ax1.scatter(atmTemp, maxOzone, s=pressure, c=starTemp, cmap='plasma',\
vmin=2000, vmax=6000)
cbar = pp.colorbar(scatter1)
cbar.set_label('Stellar Temperature [K]')
pp.legend(*scatter1.legend_elements("sizes", num=None), title='$P/P_{Earth}$ [atm]',\
 loc='upper left')
pp.xlabel('Atmosphere Temperature [K]')
ax1.set_ylabel('Maximum Ozone Abundance\n[mol/cm$^3$]')
pp.savefig('abundanceMax.png')

fig, ax1 = pp.subplots(1, sharex=True)
scatter1 = ax1.scatter(atmTemp, altMax, s=pressure, c=starTemp, cmap='plasma',\
 vmin=2000, vmax=6000)
cbar = pp.colorbar(scatter1)
cbar.set_label('Stellar Temperature [K]')
pp.legend(*scatter1.legend_elements("sizes", num=None), title='$P/P_{Earth}$ [atm]',\
 loc='upper left')
pp.xlabel('Atmosphere Temperature [K]')
ax1.set_ylabel('Altitude of Maximum Ozone Abundance\n[km]')
pp.savefig('altitudeMax.png')

fig, ax1 = pp.subplots(1, sharex=True)
scatter1 = ax1.scatter(atmTemp, totalColumn, s=pressure, c=starTemp, cmap='plasma',\
 vmin=2000, vmax=6000)
cbar = pp.colorbar(scatter1)
cbar.set_label('Stellar Temperature [K]')
pp.legend(*scatter1.legend_elements("sizes", num=None), title='$P/P_{Earth}$ [atm]',\
 loc='upper left')
pp.xlabel('Atmosphere Temperature [K]')
ax1.set_ylabel('Total Column\n[mol/cm$^3$]')
pp.savefig('totalColumn.png')

fig, (ax1, ax2, ax3) = pp.subplots(3, sharex=True)
scatter1 = ax1.scatter(atmTemp, maxOzone, s=pressure, c=starTemp, cmap='plasma',\
 vmin=2000, vmax=6000)
scatter2 = ax2.scatter(atmTemp, altMax, s=pressure, c=starTemp, cmap='plasma',\
 vmin=2000, vmax=6000)
scatter3 = ax3.scatter(atmTemp, totalColumn, s=pressure, c=starTemp, cmap='plasma',\
 vmin=2000, vmax=6000)
ax1.set_ylabel('Maximum Ozone Abundance\n [mol/cm$^3$]', fontsize=6)
pp.colorbar(scatter1)
ax2.set_ylabel('Altitude of Maximum Ozone\n[km]', fontsize=6)
ax3.set_ylabel('Total Column\n[mol/cm$^3$]', fontsize=6)
pp.xlabel('Atmosphere Temperature [K]')
pp.savefig('alldata.png')

totalInputs = len(pressure)

with open('sortedData.txt', 'w') as f:
    f.write("#Pressure, Star Temperature, Average Atm. Temp [K], Ozone Abundance [molecules/cm3], Altitude of Max [km], Total Column [molecules/cm3]\n")
    for element in range(0, totalInputs):
        f.write("{:04.1f}, {:05d}, {:05.1f}, {:13.0f}, {:05.2f}, {:14.0f}\n".format(pressure[element], starTemp[element], atmTemp[element], maxOzone[element], altMax[element], totalColumn[element]))
