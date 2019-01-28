from scipy import interpolate
import numpy as np
import settings as s
import pdb

def log_interp1d(xx, yy, kind='linear'):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, kind=kind)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

def getNOComposition(altitude):
    #NO2 and OH profiles taken from Cunnold et al. 1975

    NO2Density = [1e5,2e8,3e9,1.2e9,1e10]
    NO2Altitudes = [58,42,26,12,5]

    #linear between 0,1; parabolic between 1-3; linear 3-4
    tempden = [0]*len(altitude)
    for i in range(len(altitude)):
        if altitude[i] >= NO2Altitudes[1]:
            #linear
            tempden[iAlt] = log_interp1d(NO2Altitudes[0:2],
                NO2Density[0:2])

        elif altitude[i] < NO2Altitudes[1] and \
            altitude[i] >= NO2Altitudes[3]:
            pass #parabolic

        elif altitude[i] < NOAltitudes[3]:
            tempden[iAlt] = log_interp1d(NO2Altitudes[3:5],
                NO2Density[3:5])

    return tempden

def initializeAtmosphere(f):
    input_data= np.loadtxt(f, unpack=True) #Radiation stuff

    P=[0*i for i in s.Temperature]
    N=[0*i for i in s.Temperature]
    #Build up background atmosphere
    P[-1]=input_data[0]*1013250 #Barye
    #1 KPa = 1e5 Barye

    N[-1]=(P[-1]/s.consts['k_B']/s.Temperature[-1])
    #assume only O2 and N2
    s.vo2=input_data[1]
    mu=(32.*s.vo2+28.*(1.-s.vo2))

    s.R_star= float(input_data[2])*s.consts['R_sun']
    s.D_pl=float(input_data[4])*s.consts['au']
    s.R_pl=float(input_data[5])*s.consts['R_earth']
    s.M_pl=float(input_data[6])*s.consts['M_earth']
    s.tstar=float(input_data[3])
    s.dtOut = float(input_data[8])*86400.
    s.g=s.consts['G']*s.M_pl/s.R_pl/s.R_pl #cm s-2
    Hsca=(s.consts['k_B']*np.mean(s.Temperature)/(mu*  \
        s.consts['protonmass']*s.g))#cm

    # valori planetari e stellari

    s.P=[P[-1]*np.exp(-z/Hsca) for z in s.Altitude]
    s.N=[(N[-1]*np.exp(-z/Hsca)) for z in s.Altitude]

    s.O2=[s.vo2*i for i in s.N]
    s.N2=[(1.-s.vo2)*i for i in s.N]
    s.O3=[1e6]*len(s.N)
    s.O=[1e9]*len(s.N)
    # pp.plot(s.O2,s.Altitude)
    # pp.savefig('plot.png')

    # s.NO = getNOComposition(altitude)

    return 0
