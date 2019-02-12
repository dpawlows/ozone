#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:46:28 2018

@author: eleonoraalei
"""
import pdb
import glob
from astropy import constants as const
from astropy import units as u
import chemistry
import photo
import initAtmosphere
import numpy as np
import settings as s
import inputs
import output
from matplotlib import pyplot as pp


path_input="input/"
files=glob.glob(path_input+"*.inp")
iError = s.init()
sza = 0
iFile = 0
for f in files:
    iFile = iFile+1

    #Initialize globals
    iError = initAtmosphere.initializeAtmosphere(f)
    iError = photo.getPhotoCrosssections()

    if not usePhotoData:
        #If irraidiance is not specified, calculate BB profile
        Ftoa = photo.getBBspectrum()

    maxdiffarr=[]
    diffarr=[]
    counter=0

    diffO3=[0*i for i in s.O3]
    diffO3=np.array(diffO3)
    dO3=[0*i for i in s.O3]
    dO3=np.array(dO3)

    print("Starting main time loop (file {})".format(iFile))

    while s.totaltime.total_seconds() < s.runTime:
        ##### Main time loop #####
        if s.totaltime.total_seconds() % s.dtprint == 0:
            #Print time stamp to screen
            iError = s.printMessage()

        #get updated irradiance values
        if inputs.usePhotoData:
            Ftoa = photo.getIrradiance()

        # ATMOSPHERIC LAYERS
        tau = [0.0]*len(s.wavelengthLow)
        for iAlt in range(0,len(s.Altitude)):
            PhotoDissRate = [0.0]*s.nPhotoSpecies
            density = [s.O[iAlt],s.O2[iAlt],s.O3[iAlt],s.NO2[iAlt],
            s.NO[iAlt]]

            for iWave in range(len(s.wavelengthLow)):
                for species in s.PhotoSpecies:
                    speciesIndex = s.PhotoSpecies.index(species)

                    if iAlt > 0:
                        densityup = [s.O[iAlt-1],s.O2[iAlt-1],s.O3[iAlt-1],s.NO2[iAlt-1],
                        s.NO[iAlt-1]]
                        #tau is zero at the top of the atmosphere
                        tau[iWave] = tau[iWave] + \
                        densityup[s.cSpeciesNames.index(species)]\
                            *s.PhotoDissociationCrosssections[\
                            speciesIndex][iWave]\
                            * (s.Altitude[iAlt-1]-s.Altitude[iAlt])

                intensity = Ftoa[iWave] * \
                    np.exp(-tau[iWave]/np.arccos(sza))

                for species in s.PhotoSpecies:
                    #Need two species loops to get tau first
                    speciesIndex = s.PhotoSpecies.index(species)

                    PhotoDissRate[speciesIndex] = \
                        PhotoDissRate[speciesIndex] + \
                            intensity *\
                            s.PhotoDissociationCrosssections[\
                            speciesIndex][iWave]

            y =\
            chemistry.calcChemistry(density,PhotoDissRate,s.N[iAlt],
                s.Temperature[iAlt],inputs.tstep.total_seconds(),\
                chemsolver=inputs.chemsolver,\
                iAlt=iAlt)

            s.O2[iAlt]=y[s.iO2]
            s.O3[iAlt]=y[s.iO3]
            s.O[iAlt]=y[s.iO]

            if s.O2[iAlt] < 0 or s.O[iAlt] < 0 or s.O3[iAlt] < 0:
                print('---- ErrorNegative density in main...')
                print('----Time: {}'.format(s.totaltime))
                print('----iStep: {}'.format(s.istep))
                print('----iAlt: {}'.format(iAlt))
                pdb.set_trace()

        iError
        s.totaltime += inputs.tstep
        s.istep += 1

        if s.totaltime.total_seconds() % inputs.dtOut <\
            inputs.tstep.total_seconds():
            iError = output.output()


    iError = s.finalize()
