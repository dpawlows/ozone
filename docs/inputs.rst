=======
Inputs
=======
Each of the inputs are listed below.  The order of each input type
**does not** matter.  EOM searches for the specific headers (including
the # sign) and reads the information accordingly.  The order of
the parameters within each type **does** need to be as written here.

#DTOUT
seconds

The amount of time in between successive outputs.  Output
files are in ascii and stored in the data directory.

#TEND
ndays

The number of days to simulate

#TSTEP
nseconds

The time step.  Currently, a time step larger than 900s should
not be used.

#USEPHOTODATA
T or F
irradiancefile

If T, uses observations of the solar irradiance at the
top of the atmosphere.  EOM then reads the irradiance
from the file specified

If F, EOM assumes blackbody and calculates the irradiance
based on the parameters set using the #RADIATIONPARAMETERS flag.

#RADIATIONPARAMETERS
Pressure/Pressure_Earth
O2MixingRatio
Rstar/Rsun
Tstar
planetDistance/1AU
rPlanet/rEarth
mPlanet/mEarth
T_effective

Specifies basics of the planet, including fundamental values used
for setting up the background atmosphere and radiation environment.
Note that several of these are normalized values (e.g. a value of
1 for mPlanet/mEarth means the simulated planet has the same mass as
Earth)

#CHEMSOLVER
simple or explicit

Simple chemical solver assumes photochemical equilibrium for certain
species.  Explicit makes no such assumption. Simple is recommended.
 If explicit is used,
it is necessary for the user to determine a time step that is
appropriate for the speed of the chemistry. 
