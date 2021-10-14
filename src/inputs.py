"""Methods required for reading user input file.
"""
from datetime import datetime,timedelta
import settings as s

iError = 0

def readint(f):
    """Read a floting point variable and return."""

    temp = f.readline()
    t = temp.split()

    try:
        return int(t[0]),0
    except:
        return None,1

def readfloat(f):
    """Read a floting point variable and return."""

    temp = f.readline()
    t = temp.split()

    try:
        return float(t[0]),0
    except:
        return None,1

def readlogical(f):
    """Read a logical and return."""

    temp = f.readline()
    if temp.strip().upper()[0] == "T":
        value = True
    else:
        value = False

    return value,0

def readchar(f):
    """Read a character string and return."""

    temp = f.readline()
    t = temp.strip().lower().split()[0]
    if t == None:
        iError = 1
    else:
        iError = 0

    return t,iError

def readInputData(file):
    """Read user input file.

    User file(s) are passed as input, but are expected to be
    located in inputs/\*inp.
    """
    global usePhotoData, pressure, O2mixingratio, CO2mixingratio, rPlanet, avgmu
    global rStar, tStar, distancePlanet, massPlanet, tempPlanet
    global dtOut, tstep, tEnd, chemsolver, sza, photoFile
    global startTime, endTime, eccentricity,nDaysInYear
    global doSteadyState, steadyCondition
    global cOutputType,nOutputTypes
    global temperatureScheme, tempScaleFactor, albedo,emission
    global latitude,longitude, tilt, nHoursPerDay, SZA

    temperatureScheme = ""
    f = open(file,'r')

    iError = 0
    for line in f:
        # Location option isn't currently used
        # if line.strip().upper() == "#LOCATION":
        #     latitude,iErr = readfloat(f)
        #     iError += iErr
        #     longitude,iErr = readfloat(f)
        #     iError += iErr

        # if iError > 0:
        #     print("Error in readInputData")
        #     print("#LOCATION")
        #     print("Float    (latitude)")
        #     print("Float    (longitude)")
        #     exit(iError)

        if line.strip().upper() == "#STAR":
            tStar,iErr = readfloat(f)
            iError += iErr
            rStar,iErr = readfloat(f)
            iError += iErr

            if iError > 0:
                print("Error in readInputData")
                print("#STAR")
                print("Float    (Temperature)")
                print("Float    (r/rSun)")
                exit(iError)

        if line.strip().upper() == "#PLANET":
            distancePlanet,iErr = readfloat(f)
            s.orbitalDistance = distancePlanet
            iError = max(iError,iErr)
            rPlanet,iErr = readfloat(f)
            iError = max(iError,iErr)
            if rPlanet < 0.001:
                iError += 1
            massPlanet,iErr = readfloat(f)
            iError = max(iError,iErr)
            if massPlanet < 0.001:
                iError += 1
            tilt,iErr = readfloat(f)
            iError = max(iError,iErr)
            eccentricity,iErr = readfloat(f)
            iError = max(iError,iErr)
            nDaysInYear,iErr = readfloat(f)
            s.nSecondsPerYear = nDaysInYear * 86400.
            iError = max(iError,iErr)
            SZA,iErrr = readfloat(f)
            iError = max(iError,iErr)

            # nHoursPerDay,iError = readfloat(f)
            # iError = max(iError,iErr)
            # s.nSecondsPerYear = nDaysInYear * nHoursPerDay * 3600.

            if iError > 0:
                print("Error in readInputData")
                print("#PLANET")
                print("Float        (distance/1AU)")
                print("Float        (r/rEarth)")
                print("Float        (m/mEarth)")
                print("float        (eccentricty)")
                print("float        (tilt, degrees)")
                print("float        (ndaysinyear)")
                print("float        (SZA)")

                exit(iError)

        if line.strip().upper() == "#TEMPERATURE":
            temperatureScheme,iError = readchar(f)

            if temperatureScheme.lower() == "scaled":
                #Simple scaling based on position of planet
                tempScaleFactor, iError = readfloat(f)

            if temperatureScheme.lower() == "isothermal":

                #isothermal...
                temperature, iError = readfloat(f)
                s.initTemperature = \
                    [temperature for i in s.initTemperature]

            if temperatureScheme.lower() == "equilibrium":
                #calculate equilibrium temperature
                albedo, iError = readfloat(f)
                emission, iError = readfloat(f)
            if temperatureScheme.lower() == "avgscaled":
                #calculate equilibrium temperature
                albedo, iError = readfloat(f)
                emission, iError = readfloat(f)
            if temperatureScheme.lower() == "input":

                pass
            if iError > 0:
                print("Error in readInputData")
                print("#TEMPERATURE")
                print("Issue with formatting.  See docs.")
                exit(iError)


        if line.strip().upper() == "#ATMOSPHERE":
            pressure,iErr = readfloat(f)
            iError = max(iError,iErr)
            O2mixingratio,iErr = readfloat(f)
            iError = max(iError,iErr)
            CO2mixingratio,iErr = readfloat(f)
            iError = max(iError,iErr)
            avgmu,iErr = readfloat(f)
            iError = max(iError,iErr)
            if iError > 0:
                print("Error in readInputData")
                print("#RADIATIONPARAMETERS")
                print("Float    (P/PEarth)")
                print("Float    (O2mixingratio)")
                print("Float    (CO2mixingratio)")
                print("Float    (avg molecular weight atmosphere)")
                exit(iError)


        if line.strip().upper() == "#OUTPUT":
            dtOut,iError = readfloat(f)
            nOutputTypes,iError = readint(f)
            cOutputType = []
            for i in range(nOutputTypes):
                temp,iError = readchar(f)
                cOutputType.append(temp)
            if iError > 0:
                print("Error in readInputData")
                print("#OUTPUT")
                print("float (seconds)")
                print("nOutputTypes")
                print("outputtype1")
                print("outputtype2")
                print("...")
                exit(iError)

        if line.strip().upper() == "#USEPHOTODATA":
            usePhotoData,iError= readlogical(f)
            if usePhotoData:
                photoFile,iError = readchar(f)

            if iError > 0:
                print("Error in readInputData")
                print("#USEPHOTODATA")
                print("logical")
                print("character (filename)")
                exit(iError)

        if line.strip().upper() == "#TSTEP":

            tstep,iError= readfloat(f)
            try:
                tstep = timedelta(seconds=tstep)
            except:
                iError += 1

            if iError > 0:
                print("Error in readInputData")
                print("#TSTEP")
                print("float (s)")
                exit(iError)


        if line.strip().upper() == "#CHEMISTRY":
            chemsolver,iError = readchar(f)
            if  chemsolver != "simple" \
                and chemsolver != "backwardeuler":
                iError = 1
                print(chemsolver)

            if iError > 0:
                print("Error in readInputData")
                print("#CHEMISTRY")
                print("simple, backwardEuler, explicit (watch the time step for explicit!)")
                exit(iError)

        if line.strip().upper() == "#TSTART":
            syear,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            smonth,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            sday,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            shour,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            sminute,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            ssecond,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)

            startTime = datetime(syear,smonth,sday,\
                shour,sminute,ssecond)
            s.runTime = startTime


        if line.strip().upper() == "#TEND":
            syear,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            smonth,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            sday,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            shour,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            sminute,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)
            ssecond,iError = readint(f)
            if iError == 1:
                print("Error in readInputData")
                print("#TSTART")
                print("year\nmonth\nday\nhour\nminute\nsecond")
                exit(iError)

            endTime = datetime(syear,smonth,sday,\
                shour,sminute,ssecond)

        if line.strip().upper() == "#STEADYSTATE":
            doSteadyState,iError= readlogical(f)
            if iError == 1:
                print("Error in readInputData")
                print("#STEADYSTATE")
                print("logical")
                exit(iError)
            if doSteadyState:
                #stop model if change is less than some percent
                steadyCondition, iError = readfloat(f)

                if iError == 1:
                    print("Error in readInputData")
                    print("#STEADYSTATE")
                    print("logical")
                    print("percentError")
                    exit(iError)

                steadyCondition = steadyCondition * 0.01

    f.close()

    return iError
