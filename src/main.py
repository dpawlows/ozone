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
j=0
for f in files:
    j=j+1

    #Specify scaling for photo crosssections.
    phiO2_200=1.
    phiO2_400=0.
    phiO3_200=1.
    phiO3_400=1.
    phiNO2_200 = 1.
    phiNO2_400 = 1

    #Initialize globals
    iError = initAtmosphere.initializeAtmosphere(f)
    iError = photo.getPhotoCrosssections()

    #Set irraidiance based on data or blackbody
    if inputs.usePhotoData:
        # flnew = s.irradiance.copy()
        # flnew = flnew[0,:]
        # # flnew = np.zeros((len(s.wavelengthLow)))
        # for iWave in range(len(s.wavelengthLow)):
        #     PhotonEnergy = 6.626e-34*2.998e8 /  \
        #     ((s.wavelengthLow[iWave]+s.wavelengthHigh[iWave])/ \
        #     2.*1.0e-9)
        #
        #     #W/m2/nm to Photons/s/cm2/nm
        #     flnew[iWave] =  flnew[iWave]/PhotonEnergy/100.0**2
        #
        # Ftoa = flnew
        pass

    else:
        wavelengths = []
        for i in range(len(s.wavelengthLow)):
            wavelengths.append(s.wavelengthLow[i])
            wavelengths.append(s.wavelengthHigh[i])

        wavelengths = wavelengths*u.nm
        pdb.set_trace()
        flux_lam = photo.blackbody_lambda(wavelengths, s.tstar)
        fl=flux_lam.to(u.W / u.nm / u.m**2  / u.sr)
        fl=fl*(s.R_star/s.D_pl)**2*u.sr*3.1415

        # # NORMALIZATION FACTOR
        fl[0:1]=fl[0:1]/10. #For Earth since BB isn't perfect

        phfl=fl*wavelengths/const.h.decompose()/const.c
        phfl=phfl.decompose()
        phfl=phfl.to(1/u.s/u.cm**2/u.nm)

        phfl200=(phfl[0]+phfl[1])/2.*(wavelengths[1]-wavelengths[0])
        phfl400=(phfl[2]+phfl[3])/2.*(wavelengths[3]-wavelengths[2])

        Ftoa200=phfl200.value
        Ftoa400=phfl400.value
        Ftoa = [Ftoa200,Ftoa400]


    maxdiffarr=[]
    diffarr=[]
    counter=0

    diffO3=[0*i for i in s.O3]
    diffO3=np.array(diffO3)
    dO3=[0*i for i in s.O3]
    dO3=np.array(dO3)

    print("Starting main time loop (file {})".format(j))

    while s.totaltime.total_seconds() < s.runTime:
        ##### Main time loop #####
        if s.totaltime.total_seconds() % s.dtprint == 0:
            #Print time stamp to screen
            iError = s.printMessage()

        #get updated irradiance values
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
