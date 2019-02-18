.. _inputs:

=======
Inputs
=======
Each of the inputs are listed below.  The order of each input type
**does not** matter.  EOM searches for the specific headers (including
the # sign) and reads the information accordingly.  The order of
the parameters within each type **does** need to be as written here.
An example input file is created upon installation.  The user
should examine this file and modify as necessary prior to
running the code.

**Required: #DTOUT**::

  #DTOUT
  seconds

The amount of time in between successive outputs.  Output
files are in ascii and stored in the data directory.

**Required: #TSTART**::

  #TSTART
  year
  month
  day
  hour
  minute
  second

The end time of the simulation

**Required: #TEND**::

  #TEND
  year
  month
  day
  hour
  minute
  second

The end time of the simulation

**Required: #TSTEP**::

  #TSTEP
  nseconds

The time step.  Currently, a time step larger than 1200s should
not be used.  However, the time step depends on the irradiance
specified in the model.  Testing should be done when varying
the irradiance to determine the appropriate time step.

.. _photodata:

**Required: #USEPHOTODATA**::

  #USEPHOTODATA
  T or F
  irradiancefile

If T, observations of the solar irradiance at the
top of the atmosphere are used.  EOM then reads the irradiance
from the file specified

If F, EOM assumes blackbody and calculates the irradiance
based on the parameters set using the #RADIATIONPARAMETERS flag
and the irradiance file is ignored.

**Required: #ATMOSPHERE**::

  #ATMOSPHERE
  Pressure/Pressure_Earth
  O2MixingRatio


Specifies fundamental values used for setting up the background
atmosphere.

**Required: #STAR**::

  #STAR
  Rstar/Rsun
  Tstar

**Required: #PLANET**::

  #PLANET
  planetDistanceAU
  rPlanet/rEarth
  mPlanet/mEarth
  eccentricity
  nDaysInYear
  T_effective

Fundamental parameters to set up the radiation environment
and determine orbital distance.  Note: if the irradiance is
specified at the top of the atmosphere, an eccentricity of
1 should be used.

**Required: #CHEMSOLVER**::

  #CHEMSOLVER
  simple or explicit

Simple chemical solver assumes photochemical equilibrium for certain
species.  Explicit makes no such assumption. Simple is recommended.
If explicit is used,
it is necessary for the user to determine a time step that is
appropriate for the speed of the chemistry.
