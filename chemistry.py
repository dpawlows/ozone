import settings as s
import numpy as np
import glob
from matplotlib import pyplot as pp
import pdb

def updateRates(temp):
    #cm3/s or cm6/s
    s.kO2_O=6e-34*(temp/(298))**(-2.3)
    #s.kO2_O=6.6e-35*np.exp(510/temp)
    s.kO3_O=8e-12*np.exp(-2061/temp)
    #s.kO3_O=1.9e-11*np.exp(-2300/temp)
    s.kCl_O3=2.8e-11*np.exp(-2100/temp)
    s.kCl_O=3e-11*np.exp(+600/temp)
    s.kBr_O3=1.7e-11*np.exp(-6600/temp)
    s.kBr_O=1.9e-11*np.exp(+1900/temp)

    #NO2 chemistry rates from Cunnold et al. 1975
    s.kNO_O3 = 1.7e-12*np.exp(-1310/temp)
    s.kNO2_O = 9.1e-12

    return 0

def func(u,PhotoDissRate,N):
    #full explicity chemistry
    debug = 0
    sources = np.zeros((s.nMajorSpecies))
    losses = np.zeros((s.nMajorSpecies))

    if debug > 0:
        print('iO: {}'.format(s.iO))
        print('iO2: {}'.format(s.iO2))
        print('iO3: {}'.format(s.iO3))
        print('iNO2: {}'.format(s.iNO2))
        print('iNO3: {}'.format(s.iNO3))
    ####### Photochemistry #####################
    ### O3 + hv -> O2 + O
    r = PhotoDissRate[s.iPhotoO3]*u[s.iO3]
    sources[s.iO2] += r
    sources[s.iO] += r
    losses[s.iO3] += r
    if debug > 0:
        print("O3 + hv -> O2 + O")
        print(sources)
        print(losses)
    ### O2 + hv -> O + O
    r = PhotoDissRate[s.iPhotoO2]*u[s.iO2]
    sources[s.iO] += 2*r
    losses[s.iO2] += r
    if debug > 0:
        print("O2 + hv -> O + O")
        print(sources)
        print(losses)
    ###NO2 + hv -> NO + O
    r = PhotoDissRate[s.iPhotoNO2] * u0[s.iNO2]
    losses[s.iNO2] += r
    sources[s.iO] += r
    source[s.iNO] += r
    if debug > 0:
        print("NO2 + hv -> NO + O")
        print(sources)
        print(losses)

    ####### Exchange #############################
    ### O+O2+M -> O3 + M
    r = s.kO2_O*u[s.iO]*u[s.iO2]*N
    losses[s.iO] += r
    losses[s.iO2] +=r
    sources[s.iO3] += r
    if debug > 0:
        print("O+O2+M -> O3 + M")
        print(sources)
        print(losses)
    ### O+O3 -> 2O2
    r = s.kO3_O*u[s.iO]*u[s.iO3]
    sources[s.iO2] += 2*r
    losses[s.iO] += r
    losses[s.iO3] += r
    if debug > 0:
        print("O+O3 -> 2O2")
        print(sources)
        print(losses)
    ######NO Chemistry
    ###NO + O3 -> NO2 + O2
    r = s.kNO_O3*u0[s.iNO2]*u0[s.iO3]
    sources[s.iNO2] += r
    losses[s.iO3] += r
    losses[s.iNO] += r
    sources[s.iO2] += r
    if debug > 0:
        print("NO+O3 -> NO2 + O2")
        print(sources)
        print(losses)
    ###NO2 + O -> NO + O2
    r = s.kNO2_O*u[s.iNO2]*u[s.iO]
    losses[s.iNO2] += r
    losses[s.iO] += r
    sources[s.iNO] += r
    sources[s.iO2] += r
    if debug > 0:
        print("NO2 + O -> NO + O2")
        print(sources)
        print(losses)
    ### Cl + O3 -> ClO + O2
    r = s.kCl_O3*u[s.iO3]*s.cl
    sources[s.iO2] += r
    losses[s.iO3] += r
    if debug > 0:
        print("Cl + O3 -> ClO + O2")
        print(sources)
        print(losses)
    ### Cl + O -> ClO
    r = s.kCl_O*s.cl*u[s.iO]
    losses[s.iO] += r
    if debug > 0:
        print("Cl + O -> ClO")
        print(sources)
        print(losses)

    ### Br + O3 -> BrO + O2
    r = s.kBr_O3*s.br*u[s.iO2]
    sources[s.iO2] += r
    losses[s.iO3] += r
    if debug > 0:
        print("Br + O3 -> BrO + O2")
        print(sources)
        print(losses)
    ### Br + O -> BrO
    r = s.kBr_O*s.br*u[s.iO]
    losses[s.iO] += r
    if debug > 0:
        print("Br + O -> BrO")
        print(sources)
        print(losses)
    ################################

    if debug > 0:
        print(sources-losses)
    return sources-losses



def calcChemistry(u0,iAlt,N,temp,dt,chemsolver='simple'):
    iError = updateRates(temp)
    #Chose chemical solver
    if chemsolver == "explicit":
        update = u0 + dt*func(u0,PhotoDissRate,N)

    elif chemsolver == "simple":
        #4 Reactions with O in equilibrium with o3
        sources = [0]*s.nMajorSpecies
        losses = [0]*s.nMajorSpecies
        update = [0]*s.nMajorSpecies

        #equilibrium between O and O3:
        update[s.iO] = (2*s.O2PhotoDissociationRate[iAlt]*u0[s.iO2] + \
        s.O3PhotoDissociationRate[iAlt]*u0[s.iO3]) \
        / (s.kO2_O*u0[s.iO2]*N+s.kO3_O*u0[s.iO3])
        # print(s.O2PhotoDissociationRate[iAlt],\
        # s.O3PhotoDissociationRate[iAlt])
        ### 03 chemistry
        ###O2 + O + M -> O3 + M
        r = s.kO2_O*u0[s.iO2]*update[s.iO]*N
        sources[s.iO3] += r
        # if s.istep == 750:
        # print(r,u0[s.iO2],update[s.iO],N)
        #losses[s.iO] += r

        ###O3 + hv -> O2 + O
        r = s.O3PhotoDissociationRate[iAlt]*u0[s.iO3]
        losses[s.iO3] += r
        #sources[s.iO] += r
        # if s.istep == 750:
        # print(r)
        ###O3 + O -> 2O2
        r = s.kO3_O*u0[s.iO3]*update[s.iO]
        losses[s.iO3] += r
        #losses[s.iO] += r
        # if s.istep == 750:
        # print(r)
        # if s.istep == 750:
        # pdb.set_trace()
        update[s.iO3] = u0[s.iO3]+dt*(sources[s.iO3] - losses[s.iO3])


        #Any change in O or O3 must have come from O2
        update[s.iO2] = u0[s.iO2]# - (update[s.iO]-u0[s.iO]) - \
        #     (update[s.iO3]-u0[s.iO3])

        # #No change for these:
        # update[s.iNO2] = u0[s.iNO2]
        # update[s.iNO] = u0[s.iNO]

    else:
        print('No chemistry solver specified.')
        print('Stopping in chemistry.py')
        exit(1)

    return update
