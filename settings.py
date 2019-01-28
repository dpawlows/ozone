import pdb
from astropy import constants as const
from astropy import units as u
import numpy as np

def init():
    global F200, F400, JO2, JO3, kO2_O, kO3_3, kCl_O3, kCl_O, kBr_O, kBr_O3
    global br, cl, iO, iO2, iO3, Temperature,nLayers, nMajorSpecies
    global O2PhotoDissociationRate, O3PhotoDissociationRate
    global cSpeciesNames, nWavelengths,wavelenghtLow,wavelengthHigh
    global g, consts,Altitude,istep

    Altitude,Temperature= np.loadtxt('input/ustspline.txt', usecols=(0,1), unpack=True)

    Altitude = Altitude[::-1]*1e5
    Temperature = Temperature[::-1]

    nMajorSpecies = 3
    nWavelenghts = 2
    wavelengthLow = [199.5e9,399.5e9]
    wavelengthHigh = [200.5e9,400.5e9]
    nLayers = len(Altitude)
    #Fluxes
    F200=[0]*nLayers
    F400=[0]*nLayers
    g = 0.0
    O2PhotoDissociationRate = [0]*nLayers
    O3PhotoDissociationRate = [0]*nLayers
    #Dissocation rates
    O2PhotoDissociationRate = [0]*nLayers
    O3PhotoDissocationRate = [0]*nLayers

    cl=3e-9
    br=2e-11

    iO = 0
    iO2 = 1
    iO3 = 2
    cSpeciesNames = ['']*nMajorSpecies
    cSpeciesNames[iO] = 'O'
    cSpeciesNames[iO2] = 'O2'
    cSpeciesNames[iO3] = 'O3'

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

def plotComposition(density,altitude,time=None):
    from matplotlib import pyplot as pp
    fig=pp.figure()
    ax = fig.add_subplot(221)
    pp.tight_layout()

    for iSpecies in range(nMajorSpecies):
        ax.semilogx(density[iSpecies,:],altitude/100000.,
            label=cSpeciesNames[iSpecies])
        pp.xlabel('Density (#/cm3)')
        pp.ylabel('Altitude (km)')
        pp.legend(loc='upper left',frameon=False)
        pp.text(.02, 1.02, '{}d'.format(time/86400.),
             transform=ax.transAxes)
        pp.savefig('plot.png')

    return 0
