"""Methods required for reading user input file.
"""

import settings as s
import pdb

iError = 0

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
    global usePhotoData, pressure, O2mixingratio, rPlanet
    global rStar, tStar, distancePlanet, massPlanet, tempPlanet
    global dtOut, tstep, tEnd, chemsolver, sza, photoFile

    f = open(file,'r')

    iError = 0
    for line in f:
        if line.strip().strip().upper() == "#RADIATIONPARAMETERS":
            pressure,iErr = readfloat(f)
            iError = max(iError,iErr)
            O2mixingratio,iErr = readfloat(f)
            iError = max(iError,iErr)
            rStar,iErr = readfloat(f)
            iError = max(iError,iErr)
            tStar,iErr = readfloat(f)
            iError = max(iError,iErr)
            distancePlanet,iErr = readfloat(f)
            iError = max(iError,iErr)
            rPlanet,iErr = readfloat(f)
            iError = max(iError,iErr)
            massPlanet,iErr = readfloat(f)
            iError = max(iError,iErr)
            tempPlanet,iErr = readfloat(f)
            iError = max(iError,iErr)

            if iError > 0:
                print("Error in readInputData")
                print("#RADIATIONPARAMETERS")
                exit(iError)

        if line.strip().upper() == "#DTOUT":
            dtOut,iError = readfloat(f)

            if iError > 0:
                print("Error in readInputData")
                print("#DTOUT")
                print("float (seconds)")
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
            if iError > 0:
                print("Error in readInputData")
                print("#TSTEP")
                print("float (s)")
                exit(iError)

        if line.strip().upper() == "#TEND":
            tEnd,iError= readfloat(f)

            if iError > 0:
                print("Error in readInputData")
                print("#TSTEP")
                print("float (days)")
                exit(iError)

        if line.strip().upper() == "#CHEMISTRY":
            chemsolver,iError = readchar(f)
            if  chemsolver != "simple" and chemsolver != "explicit":
                iError = 1

            if iError > 0:
                print("Error in readInputData")
                print("#CHEMISTRY")
                print("simple or explicit (watch the time step for explicit!)")
                exit(iError)


    f.close()

    return iError
