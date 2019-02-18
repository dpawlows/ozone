import settings as s
import inputs
from numpy import cos, pi
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
    theta = (theta - s.longitudeOfPerihelion)
    if theta < 0:
        theta += 360

    distance = inputs.distancePlanet*(1-inputs.eccentricity**2)/ \
        (1+inputs.eccentricity*cos(theta*pi/180.))
    return distance
