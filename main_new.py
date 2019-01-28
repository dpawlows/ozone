import pdb

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:46:28 2018

@author: eleonoraalei
"""



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


def func(u,x):
    varo=(-k1[i1]*u[0]*u[1]*N[i1]-k2[i1]*u[0]*u[2]+JO3[i1]*u[2]+2*JO2[i1]*u[1]-k4[i1]*cl*u[0]-k6[i1]*br*u[0])
    varo2=(-k1[i1]*u[0]*u[1]*N[i1]+2*k2[i1]*u[0]*u[2]+JO3[i1]*u[2]-JO2[i1]*u[1]+k3[i1]*cl*u[2]+k4[i1]*cl*u[0]+k5[i1]*br*u[2]+k6[i1]*br*u[0])
    varo3=(k1[i1]*u[0]*u[1]*N[i1]-k2[i1]*u[0]*u[2]-JO3[i1]*u[2]-k3[i1]*cl*u[2]-k5[i1]*br*u[2])

    return [varo,varo2,varo3]

#%%

import glob
import numpy as np
from scipy.integrate import odeint
from astropy import units as u
from astropy import constants as const

R_sun=const.R_sun.to(u.cm).value #cm
R_earth=const.R_earth.to(u.cm).value #cm
M_sun=const.M_sun.to(u.g).value #g
M_earth=const.M_earth.to(u.g).value #g
au=u.au.to(u.cm) #m
G=const.G.cgs.value# cm3 g-1 s-2
k_B=const.k_B.cgs.value #erg K-1
m_p=const.m_p.cgs.value#g
cl=3e-9
br=2e-11

path_input="./input/"


files=glob.glob(path_input+"*.inp")
j=0
for f in files:
    j=j+1
    input_data= np.loadtxt(f, unpack=True)


    I0_200=9.93e12
    I0_400=1.25e14
    '''costantS'''
    sO2_200=2.51e-23 #cm2
    sO2_400=1.27e-26 #cm2
    sO3_200=2e-18 #cm2
    sO3_400=1.37e-23 #cm2
    phiO2_200=1.
    phiO2_400=0.
    phiO3_200=1.
    phiO3_400=1.
    costheta=1.


    # valori planetari e stellari
    R_star= float(input_data[2])*R_sun
    D_pl=float(input_data[4])*au
    R_pl=float(input_data[5])*R_earth
    M_pl=float(input_data[6])*M_earth
    g=G*M_pl/R_pl/R_pl #cm s-2
    tstar=float(input_data[3])


    ''' ATMOSPHERE'''
    Z,T= np.loadtxt('input/ustspline.txt', usecols=(0,1), unpack=True)
    Z=Z[::-1]
    T=T[::-1]
    Z=Z*1e5
    P=[0*i for i in T]
    N=[0*i for i in T]
    t=len(T)-1
    P[t]=input_data[0]*1013250 #Barye
    N[t]=(P[t]/k_B/T[t])
    mu=(32.*input_data[1]+28.*(1.-input_data[1]))
    Hsca=(k_B*np.mean(T)/mu/m_p/g)#cm
    P=[P[t]*np.exp(-z/Hsca) for z in Z]
    N=[(N[t]*np.exp(-z/Hsca)) for z in Z]
    vo2=input_data[1]
    O2=[vo2*i for i in N]
    N2=[(1.-vo2)*i for i in N]
    O3=[1e6]*len(N)
    O=[1e6]*len(N)
    NO=[1e6]*len(N)


    wavelengths=[199.5,200.5,399.5,400.5]*u.nm

    flux_lam = blackbody_lambda(wavelengths, tstar)
    fl=flux_lam.to(u.W / u.nm / u.m**2  / u.sr)
    fl=fl*(R_star/D_pl)**2*u.sr*3.1415

    # NORMALIZATION FACTOR (DISABLED)
    fl[0:1]=fl[0:1]/10.

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

    F200=[0 for i in T]
    F400=[0 for i in T]
    JO2_200=[0 for i in T]
    JO2_400=[0 for i in T]
    JO3_200=[0 for i in T]
    JO3_400=[0 for i in T]
    JO2=[0 for i in T]
    JO3=[0 for i in T]
    k1=[6e-34*(i/(298))**(-2.3)for i in T]
    k2=[8e-12*np.exp(-2061/i)for i in T]
    k3=[2.8e-11*np.exp(-2100/i) for i in T]
    k4=[3e-11*np.exp(+600/i) for i in T]
    k5=[1.7e-11*np.exp(-6600/i) for i in T]
    k6=[1.9e-11*np.exp(+1900/i) for i in T]
    F200[0]=Ftoa200
    F400[0]=Ftoa400

    JO2_200[0]=sO2_200*F200[0]*phiO2_200
    JO3_200[0]=sO3_200*F200[0]*phiO3_200

    JO2_400[0]=sO2_400*F400[0]*phiO2_400
    JO3_400[0]=sO3_400*F400[0]*phiO3_400

    JO2[0]=JO2_200[0]+JO2_400[0]
    JO3[0]=JO3_200[0]+JO3_400[0]

    tstep=1
    exp=np.arange(-3.,50.,0.1)
    time=0
    diffO3=[0*i for i in O3]
    diffO3=np.array(diffO3)
    dO3=[0*i for i in O3]
    dO3=np.array(dO3)
    for j in np.arange(0,len(exp)):
        O3old=O3.copy()
        O2old=O2.copy()
        Oold=O.copy()
        tstep=3600.*10.**exp[j]
        time=time+tstep
        i1=0

        u0 = [O[0], O2[0],O3[0]]
        tempo=[0,1]
        y = odeint(func, y0=u0,t=tempo)

        o=y[-1]

        O2[i1]=(o[1]) if o[1]>0 else 0
        O3[i1]=(o[2])  if o[2]>0 else 0
        O[i1]=(o[0]) if o[0]>0 else 0

        # ATMOSPHERIC LAYERS
        for i1 in range(1,len(T)):

            tau_200=(O2[i1-1]*sO2_200+O3[i1-1]*sO3_200)*(Z[i1-1]-Z[i1])
            tau_400=(O2[i1-1]*sO2_400+O3[i1-1]*sO3_400)*(Z[i1-1]-Z[i1])

            F200[i1]=(F200[i1-1]*np.exp(-tau_200/costheta))
            F400[i1]=(F400[i1-1]*np.exp(-tau_400/costheta))
            print(Z[i1]/1e5,tau_200)

            #		ABSORPTION UPPER LAYERS
            JO2_200[i1]=(sO2_200*F200[i1]*phiO2_200)
            JO3_200[i1]=(sO3_200*F200[i1]*phiO3_200)

            JO2_400[i1]=(sO2_400*F400[i1]*phiO2_400)
            JO3_400[i1]=(sO3_400*F400[i1]*phiO3_400)


            O2PhotoDissociationRate =(JO2_200[i1]+JO2_400[i1])
            O3PhotoDissocationRate =(JO3_200[i1]+JO3_400[i1])

            u0 = [O[i1], O2[i1],O3[i1]]

            y = odeint(func, y0=u0,t=tempo)
            # pdb.set_trace()
            o=y[-1]

            O2[i1]=(o[1]) if o[1]>0 else 0
            O3[i1]=(o[2])  if o[2]>0 else 0
            O[i1]=(o[0]) if o[0]>0 else 0
            if i1 >= len(T)-10:
                pdb.set_trace()
        pdb.set_trace()
        for i in range(0,len(O3)-1):
    #        diffO3[i]=(O3[i]-O3old[i])/O3[i]
            dO3[i]=(O3[i]-O3old[i])#/O3old[i]

        if len(diffarr)<2:
            diffarr.append(1000)
            conv=10
            prova=[]
        else:
            diffarr.append(np.min(dO3) if np.max(abs(dO3)) == abs(np.min(dO3)) else np.max(dO3))
            conv=np.max(O3)/np.max(O3old)-1
            prova.append(diffarr[-1]*diffarr[-2]/abs(diffarr[-1]*diffarr[-2]))





        '''PRINT OUTPUTS'''
        print(j,conv,time/3600./365.,"%E" % max(O3),diffarr[-1])

        if abs(conv)<1E-5:
            filename="final.txt"
            print(time/tstep,'%.4E' % max(O3),diffarr[-1],conv,prova[-1])
            file_out = open(filename, 'w')
            file_out.write('Z(cm)\tT(K)\tP(Ba)\tN(cm-3)\tO2(cm-3)\tO3(cm-3)\tO(cm-3)\tN2(cm-3)\tJO2(s-1)\tJO3(s-1)\n')

            for z,t,p,n,o2,o3,o,m,jo2,jo3 in zip(Z,T,P,N,O2,O3,O,N2,JO2,JO3):
                        file_out.write(str(z) + '\t' + str(t) + '\t'+str(p/10) + '\t'+ str(n)+'\t'+str(o2) + '\t'+str(o3) + '\t'+str(o) + '\t'+str(m) +"\t"+ str(jo2) + '\t'+str(jo3) +  '\n')
            file_out.close()
            break

        if ((j+1) %50 ==0):
            filename="it"+str(j)+".txt"
            file_out = open(filename, 'w')
            file_out.write('Z(cm)\tT(K)\tP(Ba)\tN(cm-3)\tO2(cm-3)\tO3(cm-3)\tO(cm-3)\tN2(cm-3)\tJO2(s-1)\tJO3(s-1)\n')

            for z,t,p,n,o2,o3,o,m,jo2,jo3 in zip(Z,T,P,N,O2,O3,O,N2,JO2,JO3):
                        file_out.write(str(z) + '\t' + str(t) + '\t'+str(p) + '\t'+ str(n)+'\t'+str(o2) + '\t'+str(o3) + '\t'+str(o) + '\t'+str(m) +"\t"+ str(jo2) + '\t'+str(jo3) +  '\n')
            file_out.close()
