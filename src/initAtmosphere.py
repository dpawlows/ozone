"""Methods required to initialize the background atmosphereself.
"""

from scipy import interpolate
import numpy as np
import settings as s
import pdb
from matplotlib import pyplot as pp
import inputs
from photo import initIrradiance
from datetime import timedelta

def log_interp1d(xx, yy, kind='linear',fill_value=''):
    """Perform log interpolation in 1D
    """

    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(xx, logy,\
        kind=kind,fill_value=fill_value)
    log_interp = lambda zz: np.power(10.0, lin_interp(zz))

    return log_interp

def getOHComposition(altitude):
    """Retrieve vertical OH profile.

    NO2 and OH profiles taken from Cunnold et al. 1975
    """

    OHDensity = [1.2e6,1.1e7]
    OHAltitudes = [i*1e5 for i in [70,41.5]]
    tempden = [0]*len(altitude)

    func1 = log_interp1d(OHAltitudes[0:2], \
        OHDensity[0:2],fill_value="extrapolate")
    tempden = func1(altitude)

    return tempden

def getNOComposition(altitude):
    """Retrieve vertical NO2 profile.

    NO2 and OH profiles taken from Cunnold et al. 1975
    """

    NO2Density = [1e5,2e8,3e9,1.45e9,1.1e9,1.25e9,1e10]
    NO2Altitudes = [i*1.e5 for i in [58,42,26,18,12,11,5]]

    #linear between 0,1; parabolic between 1-3; linear 3-4
    tempden = [0]*len(altitude)
    funchigh = log_interp1d(NO2Altitudes[0:2],
        NO2Density[0:2],fill_value="extrapolate")
    funclow = log_interp1d(NO2Altitudes[5:7],
        NO2Density[5:7],fill_value="extrapolate")
    funcmiddle1 = \
        log_interp1d(NO2Altitudes[1:4],\
            NO2Density[1:4],kind='quadratic')
    funcmiddle2 = \
        log_interp1d(NO2Altitudes[3:6],\
            NO2Density[3:6],kind='quadratic')



    for iAlt in range(len(altitude)):

        if altitude[iAlt] >= NO2Altitudes[1]:
            tempden[iAlt] = funchigh(altitude[iAlt])

        elif altitude[iAlt] < NO2Altitudes[1] and \
            altitude[iAlt] >= NO2Altitudes[3]:
            tempden[iAlt] = float(funcmiddle1(altitude[iAlt]))

        elif altitude[iAlt] < NO2Altitudes[3] and \
            altitude[iAlt] >= NO2Altitudes[5]:
            tempden[iAlt] = float(funcmiddle2(altitude[iAlt]))

        elif altitude[iAlt] < NO2Altitudes[5]:
            tempden[iAlt] = funclow(altitude[iAlt])


    return tempden

def initializeAtmosphere(f):
    """Initialize the background atmosphere including
    create fundamental parameters like pressure, densities, temperature,
    etc.
    """

    input_data = inputs.readInputData(f)

    Pressure=[0*i for i in s.initTemperature]
    nDensity=[0*i for i in s.initTemperature]
    #Build up background atmosphere
    Pressure[-1]=inputs.pressure*1013250 #Barye
    #1 KPa = 1e5 Barye

    nDensity[-1]=(Pressure[-1]/s.consts['k_B']/\
        s.initTemperature[-1])
    #assume only O2 and N2

    mu=(32.*inputs.O2mixingratio+28.*(1.-inputs.O2mixingratio))
    s.R_star= inputs.rStar*s.consts['R_sun']
    s.D_pl=inputs.distancePlanet*s.consts['au']
    s.R_pl=inputs.rPlanet*s.consts['R_earth']
    s.M_pl=inputs.massPlanet*s.consts['M_earth']
    s.tstar=inputs.tStar
    s.g=s.consts['G']*s.M_pl/s.R_pl/s.R_pl #cm s-2

    Hsca=(s.consts['k_B']*np.mean(s.initTemperature)/(mu*  \
        s.consts['protonmass']*s.g))#cm

    # valori planetari e stellari

    s.P=[Pressure[-1]*np.exp(-z/Hsca) for z in s.Altitude]
    s.N=[(nDensity[-1]*np.exp(-z/Hsca)) for z in s.Altitude]

    s.O2=[inputs.O2mixingratio*i for i in s.N]
    s.N2=[(1.-inputs.O2mixingratio)*i for i in s.N]
    s.O3=[1e6]*len(s.Altitude)
    s.O=[1e9]*len(s.Altitude)
    s.NO2 = getNOComposition(s.Altitude)
    s.NO =[1e9]*len(s.Altitude)
    s.OH = getOHComposition(s.Altitude)
    s.totaltime=timedelta(seconds=0)
    # s.density = np.zeros((s.nLayers,s.nMajorSpecies))
    s.density = np.array([s.O,s.O2,s.O3])

    if inputs.usePhotoData:
        s.irradiance = initIrradiance()

    return 0

def checkDone(olddensity):
    if inputs.doSteadyState:

        maxdiff = np.max((s.density - olddensity)/olddensity)
        s.difflist.pop(0)
        s.difflist.append(maxdiff)

        if maxdiff < inputs.steadyCondition:
                done = True
        else:
            done = False
            if s.difflist[-1] > np.mean(s.difflist) and \
                s.istep > len(s.difflist):
                print("Simulation may be diverging?")
                print("Latest percent diff is greater than \
                the mean of the last {} times.".format(len(maxdiff)))
                print("Stopping...")
                exit(1)
    else:
        if s.runTime > inputs.endTime:
            done = True
        else:
            done = False

    return done
