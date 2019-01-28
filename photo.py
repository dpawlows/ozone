import settings as s
import numpy as np
from matplotlib import pyplot as pp
from glob import glob
import pdb

def plotCrosssections(wave,sigma,species=None):
    fig=pp.figure()
    ax = fig.add_subplot(221)
    pp.tight_layout()
    ax.loglog(wave,sigma,label=species)
    pp.xlabel('Wavelength (A)')
    pp.ylabel('Cross Section (cm2)')
    if species != None:
        pp.legend(loc='upper left')
    pp.savefig('crossection.png')
    print('Quitting after cross section plot!')
    exit(1)

def getAverageSigmas(wave,sigmas):
    '''Average the cross sections across the bins specified in
    wl and wh'''
    crosssections = []
    wl = s.wavelengthLow
    wh = s.wavelengthHigh
    if len(wl) != len(wh):
        print('----Error: New bins are not compatible in getAverageSigmas')
        exit(1)


    for i in range(len(wl)):

        try:
            indexl = next(data[0] for data in enumerate(wave)\
                if data[1] >= wl[i] and data[1] < wh[i])
        except:
            indexl = None

        try:
            indexh = next(data[0] for data in enumerate(wave)\
                if data[1] >= wh[i])
        except:
            indexh = None
        #We have the indices, now fill the crosssections
        if indexl == None:
            #no available data
            crosssections.append(0)
        elif indexh == None:
            #Data available for part of the specified spectral range
            crosssections.append(np.mean(sigmas[indexl:]))
        elif indexh-1 == indexl:
            crosssections.append(sigmas[indexl])
        else:
            crosssections.append(np.mean(sigmas[indexl:indexh-1]))
            #Data available beyond the range


    return crosssections

def getPhotoCrosssections():
    iError = 0

    #In order to automate this, the file have to be named
    #using the following convention (see input folder for example):
    #nSpecies_crosssecitons.dat
    availableCrosssectionFiles = glob("input/*crosssections.dat")

    ifile = 0
    for tSpecies in s.PhotoSpecies:
        tlist = [file for file in availableCrosssectionFiles if\
         "n"+tSpecies in file]
        if len(tlist) == 0:
            iError = 1
            message = "Missing cross section file: ",tSpecies,"\n"\
            "stopping in chemistry.py"
            return iError,message

        file = tlist[0]

        f = open(file,'r')
        started = False
        wave = []
        sigmas = []
        dissociationMap = []
        for line in f:
            if started:
                temp = line.split()
                wave.append(float(temp[0])/10.) #convert to nm
                sigmasKept = [float(temp[ireact]) for ireact in \
                    dissociationMap]
                sigmas.append(sum(sigmasKept))
            if line[0:7].lstrip() == "Lambda" and not started:
                started = True
                temp = line.split()
                for i in range(len(temp)):
                    if "/" in temp[i]:
                        dissociationMap.append(i)
        makeplot = False
        if ifile == 2 and makeplot:
            plotCrosssections(wave,sigmas,'O3')

        s.PhotoDissociationCrosssections[s.PhotoSpecies.index(tSpecies)]\
            = getAverageSigmas(wave,sigmas)
        ifile += 1

    if iError == 0:
        message = ""
    return iError,message

def getIrradiance():
    #Integrate the irradiance to get the bins that we want
    #note that the file gives values at 1nm resolution,
    #and units are .../nm, so the sum is the same as the integral

    irradianceFile = "input/sorce_ssi_l3.csv"
    wave,tempirradiance = \
        np.loadtxt(irradianceFile, usecols=(1,2),\
        delimiter=",",unpack=True,skiprows=1)

    wl = s.wavelengthLow
    wh = s.wavelengthHigh


    if len(wl) != len(wh):
        print('----Error: New bins are not compatible in getAverageSigmas')
        exit(1)

    irradiance = []

    for i in range(len(wl)):
        try:
            indexl = next(data[0] for data in enumerate(wave)\
                if data[1] >= wl[i] and data[1] < wh[i])
        except:
            indexl = None

        try:
            indexh = next(data[0] for data in enumerate(wave)\
                if data[1] >= wh[i])
        except:
            indexh = None

        #We have the indices, now fill the irradiancebins
        ###We are cheating here and not integrating since the integral
        ###is equal to the sum!!
        if indexl == None:
            #no available data
            irradiance.append(0)
        elif indexh == None:
            #Data available for part of the specified spectral range
            irradiance.append(np.sum(tempirradiance[indexl:]))

        elif indexh-1 == indexl:
            irradiance.append(tempirradiance[indexl])

        else:
            irradiance.append(np.sum(tempirradiance[indexl:indexh-1]))
            #Data available beyond the range


    return irradiance
