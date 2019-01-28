#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:46:28 2018

@author: eleonoraalei
"""
import pdb


#%% FUNCTIONS
def blackbody_nu(in_x, temperature):
    FNU = u.erg / (u.cm**2 * u.s * u.Hz)
    FLAM = u.erg / (u.cm**2 * u.s * u.AA)
    """Calculate blackbody flux per steradian, :math:`B_{\\nu}(T)`.

    .. note::

        Use `numpy.errstate` to suppress Numpy warnings, if desired.

    .. warning::

        Output values might contain ``nan`` and ``inf``.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Hz.

    temperature : number, array-like, or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} Hz^{-1} sr^{-1}`.

    Raises
    ------
    ValueError
        Invalid temperature.

    ZeroDivisionError
        Wavelength is zero (when converting to frequency).

    """
    # Convert to units for calculations, also force double precision
    with u.add_enabled_equivalencies(u.spectral() + u.temperature()):
        freq = u.Quantity(in_x, u.Hz, dtype=np.float64)
        temp = u.Quantity(temperature, u.K, dtype=np.float64)


    log_boltz = const.h * freq / (const.k_B * temp)
    boltzm1 = np.expm1(log_boltz)

    # Calculate blackbody flux
    bb_nu = (2.0 * const.h * freq ** 3 / (const.c ** 2 * boltzm1))
    flux = bb_nu.to(FNU, u.spectral_density(freq))

    return flux / u.sr  # Add per steradian to output flux unit



def blackbody_lambda(in_x, temperature):
    FNU = u.erg / (u.cm**2 * u.s * u.Hz)
    FLAM = u.erg / (u.cm**2 * u.s * u.AA)
    """Like :func:`blackbody_nu` but for :math:`B_{\\lambda}(T)`.

    Parameters
    ----------
    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Angstrom.

    temperature : number, array-like, or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} \\mathring{A}^{-1} sr^{-1}`.

    """
    if getattr(in_x, 'unit', None) is None:
        in_x = u.Quantity(in_x, u.AA)

    bb_nu = blackbody_nu(in_x, temperature) * u.sr  # Remove sr for conversion
    flux = bb_nu.to(FLAM, u.spectral_density(in_x))

    return flux / u.sr  # Add per steradian to output flux unit


# def func(u,x):
#     varo=(-k1[i1]*u[0]*u[1]*N[i1]-k2[i1]*u[0]*u[2]+JO3[i1]*u[2]+2*JO2[i1]*u[1]-k4[i1]*cl*u[0]-k6[i1]*br*u[0])
#     varo2=(-k1[i1]*u[0]*u[1]*N[i1]+2*k2[i1]*u[0]*u[2]+JO3[i1]*u[2]-JO2[i1]*u[1]+k3[i1]*cl*u[2]+k4[i1]*cl*u[0]+k5[i1]*br*u[2]+k6[i1]*br*u[0])
#     varo3=(k1[i1]*u[0]*u[1]*N[i1]-k2[i1]*u[0]*u[2]-JO3[i1]*u[2]-k3[i1]*cl*u[2]-k5[i1]*br*u[2])
#
#     return [varo,varo2,varo3]

#%%

import glob
# from scipy.integrate import odeint
from astropy import constants as const
from astropy import units as u
from chemistry import calcChemistry
import initAtmosphere
import numpy as np
import settings as s
from time import time


path_input="input/"
files=glob.glob(path_input+"*.inp")
j=0


for f in files:
    j=j+1

    #? P/Patm %O2 Rstar/Rsun Tstar Dplanet(AU) Rplanet/Rearth Mplanet/Mearth ?NotUsed

    I0_200=9.93e12
    I0_400=1.25e14
    '''costantS'''
    sO2_200=2.51e-23 #cm2
    sO2_400=1.27e-26 #cm2
    sO3_200=2e-18 #cm2
    sO3_400=1.37e-23 #cm2
    phiO2_200=1.
    phiO2_400=1.
    phiO3_200=1.
    phiO3_400=1.
    costheta=1.

    ''' ATMOSPHERE'''
 #altitude
     #temperature

    #Initialize globals
    iError = s.init()
    iError = initAtmosphere.initializeAtmosphere(f)

    wavelengths=[199.5,200.5,399.5,400.5]*u.nm

    flux_lam = blackbody_lambda(wavelengths, s.tstar)
    fl=flux_lam.to(u.W / u.nm / u.m**2  / u.sr)
    fl=fl*(s.R_star/s.D_pl)**2*u.sr*3.1415
    ####fl seems large in the 199.5 and 200.5 bin
    ####by an order of magnitude compared to data on LISIRD

    # NORMALIZATION FACTOR (DISABLED)
    fl[0:1]=fl[0:1]/10.
    print(fl)
    exit(1)
    #### VERIFICA COSTANTE SOLARE OK
    #print(np.trapz(fl,x=wavelengths))
    #result=result.to(u.W/ u.m**2)

    phfl=fl*wavelengths/const.h.decompose()/const.c
    phfl=phfl.decompose()
    phfl=phfl.to(1/u.s/u.cm**2/u.nm)


    phfl200=(phfl[0]+phfl[1])/2.*(wavelengths[1]-wavelengths[0])
    phfl400=(phfl[2]+phfl[3])/2.*(wavelengths[3]-wavelengths[2])

    Ftoa200=phfl200.value
    Ftoa400=phfl400.value

    maxdiffarr=[]
    diffarr=[]
    counter=0

    '''TOP OF ATMOSPHERE'''


    #Dissocation rates?
    JO2_200=[0 for i in s.Temperature]
    JO2_400=[0 for i in s.Temperature]
    JO3_200=[0 for i in s.Temperature]
    JO3_400=[0 for i in s.Temperature]


    #rate constants

    s.F200[0]=Ftoa200
    s.F400[0]=Ftoa400

    JO2_200[0]=sO2_200*s.F200[0]*phiO2_200
    JO3_200[0]=sO3_200*s.F200[0]*phiO3_200

    JO2_400[0]=sO2_400*s.F400[0]*phiO2_400
    JO3_400[0]=sO3_400*s.F400[0]*phiO3_400

    tstep=3600
    # exp=np.arange(-3.,50.,0.1)
    totaltime=0
    diffO3=[0*i for i in s.O3]
    diffO3=np.array(diffO3)
    dO3=[0*i for i in s.O3]
    dO3=np.array(dO3)
    s.istep = 0
    print("Starting main time loop (file {})".format(j))
    startTime = time()
    s.dtOut=100*86400.
    while totaltime < s.dtOut:

        if s.istep % 500 == 0:
            elapsedTime = time() - startTime
            print("s.istep: {}; run time: {}s; elapsed time: {:03.1f}s".format(s.istep,totaltime,elapsedTime))
        # O3old=O3.copy()
        # O2old=O2.copy()
        # Oold=O.copy()
        # tstep=3600.*10.**exp[j]


        # u0 = [O[0], O2[0],O3[0]]
        # tempo=[0,1]

        # y = odeint(func, y0=u0,t=tempo,args=(i1,N[i1]))
        # y = updateComposition(u0,0,N[0],1)
        # pdb.set_trace()
        # # pdb.set_trace()
        # o=y[-1]
        #
        # O2[i1]=(o[s.iO2]) if o[1]>0 else 0
        # O3[i1]=(o[s.iO3])  if o[2]>0 else 0
        # O[i1]=(o[s.iO]) if o[0]>0 else 0

        # ATMOSPHERIC LAYERS
        tau_200 = 0
        tau_400 = 0
        for iAlt in range(0,len(s.Altitude)):

            if iAlt > 0:
                tau_200=tau_200+(s.O2[iAlt-1]*sO2_200+s.O3[iAlt-1]*sO3_200)*(s.Altitude[iAlt-1]-s.Altitude[iAlt])
                tau_400=tau_400+(s.O2[iAlt-1]*sO2_400+s.O3[iAlt-1]*sO3_400)*(s.Altitude[iAlt-1]-s.Altitude[iAlt])

            s.F200[iAlt]=(s.F200[0]*np.exp(-tau_200/costheta))
            s.F400[iAlt]=(s.F400[0]*np.exp(-tau_400/costheta))

            #		ABSORPTION UPPER LAYERS
            JO2_200[iAlt]=(sO2_200*s.F200[iAlt]*phiO2_200)
            JO3_200[iAlt]=(sO3_200*s.F200[iAlt]*phiO3_200)

            JO2_400[iAlt]=(sO2_400*s.F400[iAlt]*phiO2_400)
            JO3_400[iAlt]=(sO3_400*s.F400[iAlt]*phiO3_400)

            s.O2PhotoDissociationRate[iAlt]=(JO2_200[iAlt]+\
                JO2_400[iAlt])
            s.O3PhotoDissociationRate[iAlt]=(JO3_200[iAlt]+\
                JO3_400[iAlt])
            # pdb.set_trace()
            u0 = [s.O[iAlt], s.O2[iAlt],s.O3[iAlt]]

            # y = odeint(func, y0=u0,t=tempo,args=(iAlt,N[iAlt]))
            y = calcChemistry(u0,iAlt,s.N[iAlt],
                s.Temperature[iAlt],tstep,chemsolver="simple")


            s.O2[iAlt]=y[s.iO2]
            s.O3[iAlt]=y[s.iO3]
            s.O[iAlt]=y[s.iO]
            print(iAlt,s.O[iAlt])
            if s.O2[iAlt] < 0 or s.O[iAlt] < 0 or s.O3[iAlt] < 0:
                print('Negative density in main...')
                print('Time: {}'.format(totaltime))
                print('iAlt: {}'.format(iAlt))
                pdb.set_trace()

            totaltime += tstep
            s.istep += 1
            # if iAlt ==0:
            # pdb.set_trace()
    density = np.array([s.O,s.O2,s.O3])
    iError = s.plotComposition(density,s.Altitude,time=totaltime)
    print("Completed in s.istep: {}; run time: {}s; elapsed time: {:03.1f}s".format(s.istep,totaltime,elapsedTime))
    print('{:g}'.format(max(s.O3)))
        # pdb.set_trace()
    #     for i in range(0,len(O3)-1):
    # #        diffO3[i]=(O3[i]-O3old[i])/O3[i]
    #         dO3[i]=(O3[i]-O3old[i])#/O3old[i]
    #
    #     if len(diffarr)<2:
    #         diffarr.append(1000)
    #         conv=10
    #         prova=[]
    #     else:
    #         diffarr.append(np.min(dO3) if np.max(abs(dO3)) == abs(np.min(dO3)) else np.max(dO3))
    #         conv=np.max(O3)/np.max(O3old)-1
    #         prova.append(diffarr[-1]*diffarr[-2]/abs(diffarr[-1]*diffarr[-2]))





        # '''PRINT OUTPUTS'''
        # print(j,conv,time/3600./365.,"%E" % max(O3),diffarr[-1])
        #
        # if abs(conv)<1E-5:
        #     filename="final.txt"
        #     print(time/tstep,'%.4E' % max(O3),diffarr[-1],conv,prova[-1])
        #     file_out = open(filename, 'w')
        #     file_out.write('Z(cm)\tT(K)\tP(Ba)\tN(cm-3)\tO2(cm-3)\tO3(cm-3)\tO(cm-3)\tN2(cm-3)\tJO2(s-1)\tJO3(s-1)\n')
        #
        #     for z,t,p,n,o2,o3,o,m,jo2,jo3 in zip(Z,T,P,N,O2,O3,O,N2,JO2,JO3):
        #                 file_out.write(str(z) + '\t' + str(t) + '\t'+str(p/10) + '\t'+ str(n)+'\t'+str(o2) + '\t'+str(o3) + '\t'+str(o) + '\t'+str(m) +"\t"+ str(jo2) + '\t'+str(jo3) +  '\n')
        #     file_out.close()
        #     break
        #
        # if ((j+1) %50 ==0):
        #     filename="it"+str(j)+".txt"
        #     file_out = open(filename, 'w')
        #     file_out.write('Z(cm)\tT(K)\tP(Ba)\tN(cm-3)\tO2(cm-3)\tO3(cm-3)\tO(cm-3)\tN2(cm-3)\tJO2(s-1)\tJO3(s-1)\n')
        #
        #     for z,t,p,n,o2,o3,o,m,jo2,jo3 in zip(Z,T,P,N,O2,O3,O,N2,JO2,JO3):
        #                 file_out.write(str(z) + '\t' + str(t) + '\t'+str(p) + '\t'+ str(n)+'\t'+str(o2) + '\t'+str(o3) + '\t'+str(o) + '\t'+str(m) +"\t"+ str(jo2) + '\t'+str(jo3) +  '\n')
        #     file_out.close()
