"""Methods required for specifying the irradiance values
at the top of the atmosphere as well as the photodissociation
cross sections.
"""

import settings as s
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as pp
from glob import glob
import inputs
import pdb
from datetime import datetime

def blackbody_nu(in_x, temperature):
    """Calculate blackbody flux per steradian, :math:`B_{\\nu}(T)`.

    .. note::

        Use `numpy.errstate` to suppress Numpy warnings, if desired.

    .. warning::

        Output values might contain ``nan`` and ``inf``.

    Parameters

    in_x : number, array-like, or `~astropy.units.Quantity`
        Frequency, wavelength, or wave number.
        If not a Quantity, it is assumed to be in Hz.

    temperature : number, array-like, or `~astropy.units.Quantity`
        Blackbody temperature.
        If not a Quantity, it is assumed to be in Kelvin.

    Returns

    flux : `~astropy.units.Quantity`
        Blackbody monochromatic flux in
        :math:`erg \\; cm^{-2} s^{-1} Hz^{-1} sr^{-1}`.


    Raises

    ValueError
        Invalid temperature.

    ZeroDivisionError
        Wavelength is zero (when converting to frequency).

    """
    FNU = u.erg / (u.cm**2 * u.s * u.Hz)
    FLAM = u.erg / (u.cm**2 * u.s * u.AA)

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
    """Like :func:`blackbody_nu` but for :math:`B_{\\lambda}(T)`.

    """
    FNU = u.erg / (u.cm**2 * u.s * u.Hz)
    FLAM = u.erg / (u.cm**2 * u.s * u.AA)

    if getattr(in_x, 'unit', None) is None:
        in_x = u.Quantity(in_x, u.AA)

    bb_nu = blackbody_nu(in_x, temperature) * u.sr  # Remove sr for conversion
    flux = bb_nu.to(FLAM, u.spectral_density(in_x))

    return flux / u.sr  # Add per steradian to output flux unit

def plotCrosssections(wave,sigma,species=None):
    '''Plotting routine for quickly checking crosssections'''
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
    settings.wavelengthHigh and settings.wavelengthLow

    wave: input array specifying the bins to average to

    sigmas: input array specifying the cross sections'''


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
    '''Fill the settings.photocrosssections arrays with
    appropriate values.

    Data taken from input/\*crosssections.dat
    files from phidrates.space.swri.edu'''
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

def initIrradiance():
    '''Integrate the irradiance to the bins specified in settings

    '''
    #Note that the file gives values at 1nm resolution,
    #and units are .../nm, so the sum is the same as the integral

    irradianceFile = inputs.photoFile
    f = open(irradianceFile,'r')

    irradiancefound = False
    s.irradianceTime = []
    tempirradiance = []
    i = 0
    for line in f:
        if irradiancefound:
                temp = line.split()
                year = int(temp[0])
                mon = int(temp[1])
                day = int(temp[2])
                hour = int(temp[3])
                min = int(temp[4])
                seco = int(temp[5])

                tempirr = [float(irr) for irr in temp[6:]]
                if len(tempirr) != nInputWavelengths:
                    print("Issue reading irradiance file")
                    print("Inconsistent number of irradiance bins")
                    print("Exiting...")
                    exit(1)

                s.irradianceTime.append(datetime(year,mon,day,hour,\
                    min,seco))
                tempirradiance.append(tempirr)

        if line[0:7] == "#HEADER":
            temp = f.readline()
            t = temp.split(":")
            try:
                nInputWavelengths = int(t[1])
                inputWavelengths = np.zeros((nInputWavelengths))
            except:
                print("Issue reading nInputWavelengths in \
                    irradiance file")
                exit(1)

        if line[0:5] == "#WAVE":
            temp = f.readline()
            t = temp.split()
            try:
                inputWavelengths = [float(wv) for wv in t]
            except:
                print("Issue reading inputWavelengths in \
                    irradiance file")
                exit(1)

        if line[0:11] == "#IRRADIANCE":
            irradiancefound = True

        i += 1
        if i > 1000 and not irradiancefound:
            print("Issue reading irradiance file")
            print("Can't locate #IRRADIANCE")
            exit(1)

    nIrradianceTimes = len(s.irradianceTime)

    wl = s.wavelengthLow
    wh = s.wavelengthHigh


    if len(wl) != len(wh):
        print('----Error: New bins are not compatible in getAverageSigmas')
        exit(1)

    irradiance = np.zeros((nIrradianceTimes,len(wl)))


    for i in range(len(wl)):
        try:
            indexl = next(data[0] for data in\
                enumerate(inputWavelengths)\
                if data[1] >= wl[i] and data[1] < wh[i])
        except:
            indexl = None

        try:
            indexh = next(data[0] for data in \
            enumerate(inputWavelengths)\
                if data[1] >= wh[i])
        except:
            indexh = None

        #We have the indices, now fill the irradiancebins
        #Integrate when necessary
        for itime in range(nIrradianceTimes):
            if indexl == None:
                #no available data
                irradiance[itime,i] = 0
            elif indexh == None:
                #Data available for part of the specified spectral
                #range
                irradiance[itime,i] =\
                    np.trapz(tempirradiance[itime][indexl:],\
                    inputWavelengths[indexl:])

            elif indexh-1 == indexl:
                irradiance[itime,i] = tempirradiance[itime][indexl]

            else:
                #Data available beyond the range
                irradiance[itime,i] =\
                    np.trapz(tempirradiance[itime][indexl:indexh-1],
                        inputWavelengths[indexl:indexh-1])


    # pdb.set_trace()
    return irradiance

def getIrradiance():
    timediff = np.array([(it-inputs.startTime).total_seconds()\
     for it in s.irradianceTime])
    f = interpolate.interp1d(timediff-timediff[0],\
        s.irradiance,axis=0)
    try:
        thisirradiance = f(s.totaltime.total_seconds())
    except ValueError:
        print("Error in getIrradiance")
        print("Values specified in irradiance file does not encompass simulation time")
        print("Exiting...")
        exit(1)
    for iWave in range(len(s.wavelengthLow)):
        PhotonEnergy = 6.626e-34*2.998e8 /  \
        ((s.wavelengthLow[iWave]+s.wavelengthHigh[iWave])/ \
        2.*1.0e-9)

        #W/m2/nm to Photons/s/cm2/nm
        thisirradiance[iWave] =  thisirradiance[iWave] \
            /PhotonEnergy/(100.0**2)


    return(thisirradiance)

def getBBspectrum():
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
