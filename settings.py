import pdb
from astropy import constants as const
from astropy import units as u
from matplotlib import pyplot as pp
import numpy as np

def init():
    global JO2, JO3, kO2_O, kO3_3, kCl_O3, kCl_O, kBr_O, kBr_O3
    global br, cl, iO, iO2, iO3, iNO2, iNO
    global Temperature,nLayers, nMajorSpecies
    global O2PhotoDissociationRate, O3PhotoDissocationRate
    global NO2PhotoDissocationRate
    global cSpeciesNames, nWavelengths,wavelengthLow,wavelengthHigh
    global g, consts, Altitude, tstar, D_pl, R_pl, R_star
    global N, O3, N2, O, P, NO
    global PhotoSpecies, PhotoDissociationCrosssections
    global iPhotoO2,iPhotoO3,iPhotoNO2, nPhotoSpecies
    global istep, dtprint, tEnd, totaltime
    global nSecondsInDay, nSecoundsInHour, nSecondsInMinute
    global density, sza

    dtprint = 3600
    istep = 0
    nSecondsInDay = 86400
    nSecoundsInHour = 3600
    nSecondsInMinute = 60

    Altitude,Temperature= np.loadtxt('input/ustspline.txt',\
     usecols=(0,1), unpack=True)

    Altitude = Altitude[::-1]*1e5
    Temperature = Temperature[::-1]

    wavelengthLow = [200,300]
    wavelengthHigh = [300,400]
    nLayers = len(Altitude)

    g = 0.0

    #Chemical species
    iO = 0
    iO2 = 1
    iO3 = 2
    iNO2 = 3
    iNO = 4
    nMajorSpecies = iNO+1
    cSpeciesNames = ['']*nMajorSpecies
    cSpeciesNames[iO] = 'O'
    cSpeciesNames[iO2] = 'O2'
    cSpeciesNames[iO3] = 'O3'
    cSpeciesNames[iNO2] = 'NO2'
    cSpeciesNames[iNO] = 'NO'


    #Dissocation Stuff
    #Different than chemical because not all majors have dissociation
    O2PhotoDissociationRate = [0]*nLayers
    O3PhotoDissocationRate = [0]*nLayers
    NO2PhotoDissocationRate = [0]*nLayers
    iPhotoO2 = 0
    iPhotoO3 = 1
    iPhotoNO2 = 2
    PhotoSpecies = ['O2','O3','NO2']
    nPhotoSpecies = len(PhotoSpecies)
    PhotoDissociationCrosssections = [0.0]*nPhotoSpecies
    cl=3e-9
    br=2e-11
    usePhotoData = False


    consts = {'R_sun':const.R_sun.to(u.cm).value, #cm
        'R_earth':const.R_earth.to(u.cm).value, #cm
        'M_sun':const.M_sun.to(u.g).value ,#g
        'M_earth':const.M_earth.to(u.g).value ,#g
        'au':u.au.to(u.cm) ,#m
        'G':const.G.cgs.value,# cm3 g-1 s-2
        'k_B':const.k_B.cgs.value, #erg K-1
        'protonmass' : const.m_p.cgs.value,#g
    }

    return 0
