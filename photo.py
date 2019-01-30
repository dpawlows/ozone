import settings as s
import numpy as np
from matplotlib import pyplot as pp
from glob import glob
import pdb

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
    
def plotCrosssections(wave,sigma,species=None):
    fig=pp.figure()
    ax = fig.add_subplot(221)
    pp.tight_layout()
    ax.loglog(wave,sigma,label=species)
    pp.xlabel('Wavelength (A)')
    pp.ylabel('Cross Section (cm2)')
    if species != None:
        pp.legend(loc='upper left')
    pp.savefig('crossection.png')
    print('Quitting after cross section plot!')
    exit(1)

def getAverageSigmas(wave,sigmas):
    '''Average the cross sections across the bins specified in
    wl and wh'''
    crosssections = []
    wl = s.wavelengthLow
    wh = s.wavelengthHigh
    if len(wl) != len(wh):
        print('----Error: New bins are not compatible in getAverageSigmas')
        exit(1)


    for i in range(len(wl)):

        try:
            indexl = next(data[0] for data in enumerate(wave)\
                if data[1] >= wl[i] and data[1] < wh[i])
        except:
            indexl = None

        try:
            indexh = next(data[0] for data in enumerate(wave)\
                if data[1] >= wh[i])
        except:
            indexh = None
        #We have the indices, now fill the crosssections
        if indexl == None:
            #no available data
            crosssections.append(0)
        elif indexh == None:
            #Data available for part of the specified spectral range
            crosssections.append(np.mean(sigmas[indexl:]))
        elif indexh-1 == indexl:
            crosssections.append(sigmas[indexl])
        else:
            crosssections.append(np.mean(sigmas[indexl:indexh-1]))
            #Data available beyond the range


    return crosssections

def getPhotoCrosssections():
    iError = 0

    #In order to automate this, the file have to be named
    #using the following convention (see input folder for example):
    #nSpecies_crosssecitons.dat
    availableCrosssectionFiles = glob("input/*crosssections.dat")

    ifile = 0
    for tSpecies in s.PhotoSpecies:
        tlist = [file for file in availableCrosssectionFiles if\
         "n"+tSpecies in file]
        if len(tlist) == 0:
            iError = 1
            message = "Missing cross section file: ",tSpecies,"\n"\
            "stopping in chemistry.py"
            return iError,message

        file = tlist[0]

        f = open(file,'r')
        started = False
        wave = []
        sigmas = []
        dissociationMap = []
        for line in f:
            if started:
                temp = line.split()
                wave.append(float(temp[0])/10.) #convert to nm
                sigmasKept = [float(temp[ireact]) for ireact in \
                    dissociationMap]
                sigmas.append(sum(sigmasKept))
            if line[0:7].lstrip() == "Lambda" and not started:
                started = True
                temp = line.split()
                for i in range(len(temp)):
                    if "/" in temp[i]:
                        dissociationMap.append(i)
        makeplot = False
        if ifile == 2 and makeplot:
            plotCrosssections(wave,sigmas,'O3')

        s.PhotoDissociationCrosssections[s.PhotoSpecies.index(tSpecies)]\
            = getAverageSigmas(wave,sigmas)
        ifile += 1

    if iError == 0:
        message = ""
    return iError,message

def getIrradiance():
    #Integrate the irradiance to get the bins that we want
    #note that the file gives values at 1nm resolution,
    #and units are .../nm, so the sum is the same as the integral

    irradianceFile = "input/sorce_ssi_l3.csv"
    wave,tempirradiance = \
        np.loadtxt(irradianceFile, usecols=(1,2),\
        delimiter=",",unpack=True,skiprows=1)

    wl = s.wavelengthLow
    wh = s.wavelengthHigh


    if len(wl) != len(wh):
        print('----Error: New bins are not compatible in getAverageSigmas')
        exit(1)

    irradiance = []

    for i in range(len(wl)):
        try:
            indexl = next(data[0] for data in enumerate(wave)\
                if data[1] >= wl[i] and data[1] < wh[i])
        except:
            indexl = None

        try:
            indexh = next(data[0] for data in enumerate(wave)\
                if data[1] >= wh[i])
        except:
            indexh = None

        #We have the indices, now fill the irradiancebins
        ###We are cheating here and not integrating since the integral
        ###is equal to the sum!!
        if indexl == None:
            #no available data
            irradiance.append(0)
        elif indexh == None:
            #Data available for part of the specified spectral range
            irradiance.append(np.sum(tempirradiance[indexl:]))

        elif indexh-1 == indexl:
            irradiance.append(tempirradiance[indexl])

        else:
            irradiance.append(np.sum(tempirradiance[indexl:indexh-1]))
            #Data available beyond the range


    return irradiance
