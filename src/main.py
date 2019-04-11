#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:46:28 2018

@author: eleonoraalei
"""
import pdb
import glob
import chemistry
import initAtmosphere
import numpy as np
import settings as s
import inputs
import output
import photo
import orbit
from matplotlib import pyplot as pp


path_input="input/"
files=glob.glob(path_input+"*.inp")
iError = s.init()
iFile = 0
for f in files:
    iFile = iFile+1

    #Initialize globals
    iErr = initAtmosphere.initializeAtmosphere(f)
    if iErr > 0:
        s.stopEOM("Error in initializeAtmosphere")
    iErr,message = photo.getPhotoCrosssections()
    if iErr > 0:
        s.stopEOM(message)

    PhotoDissRate_Alt = np.zeros((s.nPhotoSpecies,s.nLayers))
    oldDensity = np.copy(s.density)

    print("Starting main time loop (file {})".format(iFile))
    done = False
    iErr = 0
    while not done:
        ##### Main time loop #####
        if (s.runTime - inputs.startTime).total_seconds() % \
            s.dtprint == 0:
            #Print time stamp to screen
            iErr += s.printMessage(s.runTime-inputs.startTime)
        #update orbit
        iErr += orbit.getOrbitalDistance()
        sza = orbit.getSZA()

        #get updated irradiance values
        Ftoa = photo.getIrradiance()

        #update Temperature
        iErr += photo.updateTemperature()
        if iErr > 0:
            s.stopEOM("Error in main.")

        # ATMOSPHERIC LAYERS
        tau = [0.0]*len(s.wavelengthLow)


        for iAlt in range(0,len(s.Altitude)):
            PhotoDissRate = [0.0]*s.nPhotoSpecies

            for iWave in range(len(s.wavelengthLow)):
                for species in s.PhotoSpecies:
                    speciesIndex = s.PhotoSpecies.index(species)

                    if iAlt > 0:
                        #tau is zero at the top of the atmosphere
                        densityup = s.density[:,iAlt-1]
                        tau[iWave] = tau[iWave] + \
                        densityup[s.cSpeciesNames.index(species)]\
                            *s.PhotoDissociationCrosssections[\
                            speciesIndex][iWave]\
                            * (s.Altitude[iAlt-1]-s.Altitude[iAlt])

                intensity = Ftoa[iWave] * \
                    np.exp(-tau[iWave]/np.cos(sza*np.pi/180.))

                for species in s.PhotoSpecies:
                    speciesIndex = s.PhotoSpecies.index(species)

                    PhotoDissRate[speciesIndex] = \
                        PhotoDissRate[speciesIndex] + \
                            intensity *\
                            s.PhotoDissociationCrosssections[\
                            speciesIndex][iWave]

            s.PhotoDissRate_Alt[:,iAlt] = PhotoDissRate
            y =\
            chemistry.calcChemistry(s.density[:,iAlt],PhotoDissRate,\
                s.N[iAlt], s.Temperature[iAlt],\
                inputs.tstep.total_seconds(),\
                chemsolver=inputs.chemsolver,\
                iAlt=iAlt)

            s.density[:,iAlt] = y

            if np.min(s.density) < 0:
                print('---- Error: Negative density in main...')
                print('----Time: {}'.format(s.runTime))
                print('----iStep: {}'.format(s.istep))
                print('----iAlt: {}'.format(iAlt))
                exit(1)

        s.runTime += inputs.tstep
        s.istep += 1

        if (s.runTime - inputs.startTime).total_seconds() %\
         inputs.dtOut < inputs.tstep.total_seconds():
            iError = output.output()

        #check if done
        done = initAtmosphere.checkDone(oldDensity)

        oldDensity = np.copy(s.density)

    iError = s.finalize(s.runTime-inputs.startTime)
