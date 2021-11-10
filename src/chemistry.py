"""Methods related to the actual chemical processes used
in EOM
"""
import settings as s
import numpy as np
import glob
from matplotlib import pyplot as pp
from scipy import linalg
import pdb
import inputs

def updateRates(temp):
    '''Update chemical reaction rates given the temperature'''
    #cm3/s or cm6/s
    s.kO2_O=6e-34*(temp/(298))**(-2.3)
    s.kO3_O=8e-12*np.exp(-2061/temp)
    s.kCl_O3=2.8e-11*np.exp(-2100/temp)
    s.kCl_O=3e-11*np.exp(+600/temp)
    s.kBr_O3=1.7e-11*np.exp(-6600/temp)
    s.kBr_O=1.9e-11*np.exp(+1900/temp)

    #NO2 chemistry rates from Cunnold et al. 1975
    s.kNO_O3 = 1.7e-12*np.exp(-1310/temp)
    s.kNO2_O = 9.1e-12

    #CO chemistry rates from Booker et al. 2014
    s.kCO_O_CO2 = 6.5e-33*np.exp(-2180/temp)
    s.kO_O2_CO2 = 1.35e-33

    return 0

def calcsourceterms(u,PhotoDissRate,N,dt,ialt):
    '''Calculate the full time-dependent chemical sources and losses
    for the current time step explicitly.

    '''

    density = np.array(u)

    debug = 0

    sources = np.zeros((s.nMajorSpecies))
    losses = np.zeros((s.nMajorSpecies))
    if debug > 0:
        print('iO: {}'.format(s.iO))
        print('iO2: {}'.format(s.iO2))
        print('iO3: {}'.format(s.iO3))
        print('iNO2: {}'.format(s.iNO2))
        print('iNO: {}'.format(s.iNO))
        print('iCO: {}'.format(s.iCO))
        print('iCO2: {}'.format(s.iCO2))

    if s.istep == 0:
        #initialize with a reasonable O profile
        density[s.iO] = (2*PhotoDissRate[s.iPhotoO2]*density[s.iO2] + \
        PhotoDissRate[s.iPhotoO3]*density[s.iO3]) \
        / (s.kO2_O*density[s.iO2]*N+s.kO3_O*density[s.iO3])

    ####### Photochemistry #####################
    ### O3 + hv -> O2 + O
    r = PhotoDissRate[s.iPhotoO3]*density[s.iO3]
    s.userdata[1,ialt] = PhotoDissRate[s.iPhotoO3]
    s.userdata[7,ialt] = r

    sources[s.iO2] += r
    sources[s.iO] += r
    losses[s.iO3] += r
    if debug > 0:
        print("O3 + hv -> O2 + O")
        print(r)
        print(sources[2],losses[2])
    ## O2 + hv -> O + O
    r = PhotoDissRate[s.iPhotoO2]*density[s.iO2]
    s.userdata[2,ialt] = PhotoDissRate[s.iPhotoO2]
    s.userdata[8,ialt] = r

    sources[s.iO] += 2*r
    losses[s.iO2] += r
    if debug > 0:
        print("O2 + hv -> O + O")
        print(r)
        print(sources[2],losses[2])

    ## CO2 + hv -> CO + O
    r = PhotoDissRate[s.iPhotoCO2]*density[s.iCO2]
    sources[s.iCO] += r
    sources[s.iO] += r
    #losses[s.iCO2] += r
    if debug > 0:
        print("CO2 + hv -> CO + O")
        print(r)
        print(sources[2],losses[2])

    ####### Exchange #############################
    ### O+O2+M -> O3 + M
    r = s.kO2_O*density[s.iO]*density[s.iO2]*N
    s.userdata[3,ialt] = s.kO2_O
    s.userdata[9,ialt] = r
    losses[s.iO] += r
    losses[s.iO2] +=r
    sources[s.iO3] += r
    if debug > 0:
        print("O+O2+M -> O3 + M")
        print(r)
        print(sources[2],losses[2])

    ### O+O3 -> 2O2
    r = s.kO3_O*density[s.iO]*density[s.iO3]
    s.userdata[4,ialt] =  s.kO3_O
    s.userdata[10,ialt] = r
    sources[s.iO2] += 2*r
    losses[s.iO] += r
    losses[s.iO3] += r

    ### CO+O+CO2 -> 2CO2
    r = s.kCO_O_CO2*density[s.iO]*density[s.iCO]*density[s.iCO2]
    losses[s.iCO] += r
    losses[s.iO] +=r
    #sources[s.iCO2] += r
    if debug > 0:
        print("CO+O+CO2 -> 2CO2")
        print(r)
        print(sources[2],losses[2])

    ### O+O2+CO2 -> O3+CO2
    r = s.kO_O2_CO2*density[s.iO]*density[s.iO2]*density[s.iCO2]
    losses[s.iO] += r
    losses[s.iO2] +=r
    sources[s.iO3] += r
    if debug > 0:
        print("CO+O+CO2 -> O3+CO2")
        print(r)
        print(sources[2],losses[2])

    ########## Catalytic #######################
    ### Cl + O3 -> ClO + O2
    r = s.kCl_O3*density[s.iO3]*s.cl
    s.userdata[5,ialt] =  s.kCl_O3
    sources[s.iO2] += r
    losses[s.iO3] += r

    ### Br + O3 -> BrO + O2
    r = s.kBr_O3*s.br*density[s.iO3]
    s.userdata[6,ialt] =  s.kBr_O3
    sources[s.iO2] += r
    losses[s.iO3] += r

    if debug > 0:
        print("O+O3 -> 2O2")
        print(r)
        print(sources[2],losses[2])

    return sources,losses

def calcJacobian(density,PhotoDissRate,N):
    '''Get jacobian maxtrix for backward euler. '''
    #Create each element that would go in the Jacobian here, then add to the code below
    rO2O_O2 = s.kO2_O*density[s.iO2]*N #O + O2 + M -> O3 + M, track O2
    rO2O_O  = s.kO2_O*density[s.iO]*N #O + O2 + M -> O3 + M, track O
    rO3O_O3 = s.kO3_O*density[s.iO3] #O3 + O -> 2 O2, track O3
    rO3O_O  = s.kO3_O*density[s.iO]  #O3 + O -> 2 O2, track O
    rOCOCO2_COCO2 = s.kCO_O_CO2*density[s.iCO]*density[s.iCO2] #O + CO + CO2 -> 2CO2, track CO/CO2
    rOCOCO2_OCO2  = s.kCO_O_CO2*density[s.iO]*density[s.iCO2] #O + CO + CO2 -> 2CO2, track O/CO2
    rOO2CO2_O2CO2 = s.kO_O2_CO2*density[s.iO2]*density[s.iCO2] #O + O2 + CO2 -> O3 + CO2, track O2/CO2
    rOO2CO2_OCO2 = s.kO_O2_CO2*density[s.iO]*density[s.iCO2] #O + O2 + CO2 -> O3 + CO2, track O/CO2
    rClO3 = s.kCl_O3*s.cl #Cl + O3 -> ClO + O2, track Cl
    rBrO3 = s.kBr_O3*s.br #Br + O3 -> BrO + O2, track Br

    ###This defines the nxn jacobian matrix where element k_i,j is
    ###defined by the derivative of the P-L matrix.  Note that the
    ###index order depends on the structure of the density arrays
    ###and the reactions must be in the right index.
    ###Taken from the Jacobian sheet in the reactoins workbook.
    ###Each left square bracket represents the start of a new row of the Jacobian

    k = np.array([
    [-rO2O_O2-rO3O_O3-rOCOCO2_COCO2,\
    2*PhotoDissRate[s.iPhotoO2]-rO2O_O,\
    PhotoDissRate[s.iPhotoO3]-rO3O_O,\
    -rOCOCO2_OCO2],
    [2*rO3O_O3-rO2O_O2,\
    -rO2O_O-PhotoDissRate[s.iPhotoO2],\
    PhotoDissRate[s.iPhotoO3]+2*rO3O_O+rClO3+rBrO3,\
    0],
    [rO2O_O2-rO3O_O3,\
    rO2O_O,\
    -PhotoDissRate[s.iPhotoO3]-rO3O_O-rClO3-rBrO3,\
    0],
    [-rOCOCO2_COCO2,\
    0,\
    0,\
    -rOCOCO2_OCO2]
    ])

    return k

def backwardEuler(density,PhotoDissRate,N,dt,iAlt):
    '''Update the density using the backward euler method.
    Currently requires source and loss terms to be calculated
    and their derivatives.'''
    density = np.array(density)

    sources,losses = \
     calcsourceterms(density,PhotoDissRate,N,dt,iAlt)
    ynew = np.copy(density)
    yold = [1.0]*len(ynew)
    tol = 0.1
    I = np.identity(s.nMajorSpecies)
    istep = 0

    while max(abs(ynew - yold)/yold) > tol:
        #print("error thing=", (ynew-yold)/yold)

        yold = np.copy(ynew)
        #print("yold =", yold)
        #print("ynew =", ynew)
        #print("density =", density)
        #breakpoint()

        sources,losses =\
         calcsourceterms(yold,PhotoDissRate,N,dt,iAlt)

        #Make the backwards step
        g = yold[0:s.nMajorSpecies]-density[0:s.nMajorSpecies]-dt*(sources - losses)
        K = calcJacobian(yold,PhotoDissRate,N)
        #print("sources", sources, "losses", losses)
        #print(sources-losses)
        #print("g", g)
        #print("density", density)
        #print("yold-density", yold-density)
        #print("k", K)

        ynew[0:s.nMajorSpecies] = yold[0:s.nMajorSpecies] - (linalg.inv(I-K*dt)).dot(g)
        print("ynew =", ynew)
        breakpoint()
        istep = istep + 1

    return ynew


def calcChemistry(u0,PhotoDissRate,N,temp,dt,chemsolver='simple',iAlt=None,usr=None):
    """Perform the chemistry update.

    The \#CHEMISTRY input flag determines the method to use,
    simple or explicit.

    Parameters
    ----------

    u0 : array
        Input density at a single altitude for all species
    PhotoDissRate : array
        Input photo dissociation rates for all photo active
        species
    N : float
        Total number density at single altitude
    temp : float
        Temperature at single altitude
    dt : float
        Time step
    chemsolver : char, optional
        Specify which method to use to solve chemistry
    iAlt : int, optional
        Mainly used for debugging

    """

    iError = updateRates(temp)

    #Chose chemical solver

    if chemsolver == "backwardeuler":
        update = backwardEuler(u0,PhotoDissRate,N,dt,iAlt)

    elif chemsolver == "simple":
        #4 Reactions with O in equilibrium with o3
        sources = [0]*s.nMajorSpecies
        losses = [0]*s.nMajorSpecies
        update = [0]*s.nMajorSpecies

        #equilibrium between O and O3:

        update[s.iO] = (2*PhotoDissRate[s.iPhotoO2]*u0[s.iO2] + \
        PhotoDissRate[s.iPhotoO3]*u0[s.iO3]) \
        / (s.kO2_O*u0[s.iO2]*N+s.kO3_O*u0[s.iO3])

        ### 03 chemistry
        ###O2 + O + M -> O3 + M
        r = s.kO2_O*u0[s.iO2]*update[s.iO]*N
        sources[s.iO3] += r

        ###O3 + hv -> O2 + O
        r = PhotoDissRate[s.iPhotoO3]*u0[s.iO3]
        losses[s.iO3] += r
        ###O3 + O -> 2O2
        r = s.kO3_O*u0[s.iO3]*update[s.iO]
        losses[s.iO3] += r
        ### Cl + O3 -> ClO + O2
        r = s.kCl_O3*u0[s.iO3]*s.cl
        # sources[s.iO2] += r
        losses[s.iO3] += r

        ### Br + O3 -> BrO + O2
        r = s.kBr_O3*s.br*u0[s.iO3]
        # sources[s.iO2] += r
        losses[s.iO3] += r

        #update the density arrays
        update[s.iO3] = u0[s.iO3]+dt*(sources[s.iO3] - losses[s.iO3])
        update[s.iO2] = u0[s.iO2] #O2 stays constant

        #No change for these:
        # update[s.iNO2] = u0[s.iNO2]
        # update[s.iNO] = u0[s.iNO]
        # pdb.set_trace()
        usr[1] = PhotoDissRate[s.iPhotoO3]
        usr[2] = PhotoDissRate[s.iPhotoO2]
        usr[3] = s.kO2_O
        usr[4] =  s.kO3_O
        usr[5] =  s.kCl_O3
        usr[6] =  s.kBr_O3
    else:
        print('----Error: No chemistry solver specified.')
        print('----Stopping in chemistry.py')
        exit(1)


    return update
