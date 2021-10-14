
altitude = []
col1 = []
col7 = []
total = []

with open('nCO2_crosssections.dat','r') as data:
    for line in data:
        filestring = line.split()
        altitude.append(float(filestring[0]))
        col1.append(float(filestring[2]))
        col7.append(float(filestring[8]))
        totalRate = float(filestring[2]) + float(filestring[8])
        total.append(totalRate)

lengthArrays = len(total)


with open('nCO2_crosssectionsTEST.dat', 'w') as f:
    f.write('0 Branching ratio for CO2      Total        2 of 7 branches\n')
    f.write(' Lambda  Total   sCO/O1D  tCO/O\n')
    for element in range(0, lengthArrays):
        f.write('  {:5.1f} {:3.2E} {:3.2E} {:3.2E}\n'.format(altitude[element], total[element], col1[element], col7[element]))
