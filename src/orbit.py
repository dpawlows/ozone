import settings as s
import inputs
from numpy import cos, pi, sin, arctan, tan, arccos
import pdb
def getOrbitalDistance():
    """Get the orbital distance based on the simulation time
    using a constant velocity approximation."""

    #amount of time since the vernal equinox
    timediff = (s.runTime - s.equinox).total_seconds() % \
        s.nSecondsPerYear

    #angle swept out since equinox
    theta = timediff*360/s.nSecondsPerYear

    if theta > 360:
        print("Error in orbit.py")
        print("Issue with orbit angle")
        exit(1)

    #Perihelion occurs near Jan 3, so this angle is set
    #as s.longitudeOfPerihelion
    s.orbitAngle = (theta - s.longitudeOfPerihelion)
    if s.orbitAngle < 0:
        s.orbitAngle += 360

    s.orbitalDistance = inputs.distancePlanet*(1-inputs.eccentricity**2)/ \
        (1+inputs.eccentricity*cos(s.orbitAngle*pi/180.))

    return 0


def calcSZA():

    time = 12*3600

    SunDeclination =  inputs.tilt*cos(360*\
        inputs.nDaysInYear*(nDaysSinceSolstice))

    # arctan(tan(inputs.tilt*pi/180.)*\
    #     sin(s.orbitAngle*pi/180.))

    localTime = (time/3600.0 + (inputs.longitude*pi/180)*inputs.nHoursPerDay / (2*pi)) % \
        inputs.nHoursPerDay
    sinDec = sin(SunDeclination)
    cosDec = cos(SunDeclination)
    sza = arccos(sinDec*sin(inputs.latitude*pi/180) + \
        cosDec * cos(inputs.latitude*pi/180) * \
        cos(pi*(localTime - inputs.nHoursPerDay/2)/(inputs.nHoursPerDay/2)))

    print(localTime,sza*180/pi)
    pdb.set_trace()

    return 0
