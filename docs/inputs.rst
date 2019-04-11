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


**Required: #OUTPUT**::

  #OUTPUT
  dtOut         (float)
  nOutputTypes  (int)
  outputtype1   (str)
  outputtype2   (str)
  ...

Specify output.  dtOutput is the number of seconds between
writing output.  nOutputTypes specifies the number of
output types to write.  Then, each outputtype must be
specified.  See output.py for complete list.

Current
output types are ALL, PHOTO, and USER.  ALL and PHOTO
are filled with specific values of interest.

USER is
meant to be more flexible.  The user can add a set of values
to an array without needing to write a bunch of code in output.py.
See user.py for instructions.


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

The time step.  The time step depends on the irradiance
specified in the model and the chemical scheme used.
 Testing should be done when changing these things to determine the appropriate time step.  

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
  Tstar
  Rstar/Rsun

Stellar parameters for radiation calculation purposes.

**Required: #PLANET**::

  #PLANET
  planetDistanceAU
  rPlanet/rEarth
  mPlanet/mEarth
  eccentricity
  nDaysInYear
  solarZenithAngle

Fundamental parameters to set up the radiation environment
and determine orbital distance.  Note: if the irradiance is
specified at the top of the atmosphere (e.g. by data),
an eccentricity of
1 should be used.

**Required: #CHEMISTRY**::

  #CHEMISTRY
  simple or backwardEuler

Simple chemical solver assumes photochemical equilibrium for certain
species.  backwardEuler makes no such assumption and the Backward
Euler implicit scheme is used.

**Optional: #TEMPERATURE

  #TEMPERATURE
  scaled/isothermal/equilibrium/
  parameter1
  parameter2
  ...

There are a few different ways to calculate/specify the temperature.
If this option is not used, then the temperature is read from a
file located in input/ustspline.txt and held constant.  This file is required
regardless of the temperature method specified as it sets up the
vertical coordinate and the initial pressure.

Scaled: the temperature is initialized using the values in
input/ustspline.txt and then scaled depending on the orbital position.
Required parameters: 1
ScaleFactor   float; Multiply initial temperature by this factor at perihelion

isothermal: Use isothermal temperature profile
Required parameters: 1
Temperature   float; the temperature of the atmosphere

equilibrium: calculate the equilibrium temperature based on stellar and
planetary parameters (stellar parameters are specified using #STAR flag, see above).
Required parameters: 2
albedo        float; surface albedo
emmissivity    float; emmissivity (0.67 for Earth)

**Optional: #STEADYSTATE

  #STEADYSTATE
  T/F               (logical)
  stoppingCriteria  (float)

If true, run the model until a steady state has been reached,
defined by the stopping criteria (as a percent difference).
Setting this to T means that any #TEND condition that is specified
will be ignored.  If this option is not specified, or is specified as
F, then #TEND is used as the stopping condition.
