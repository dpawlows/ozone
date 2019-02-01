===============
Getting Started
===============
To install the code, clone from github an execute::

  >>install.sh

at the prompt.  This will create a few files, including
eom.py.

========
Usage
========

Running the code is as simple as invoking python to run eom.py.
However, it is probably necessary to modify the input files to suit
the needs of the user. These files are located in the **input**
directory and are titled:

input???.inp

EOM is capable of reading multiple input files to complete an
ensemble of simulations.

See :ref:`inputs` for specific input options and descriptions.

Additionally, the user may wish to alter the irradiance values that
are used.  These values need to be specified
in a single csv file
(see the :ref:`#PHOTODATA <photodata>` input)
in which the first two columns are wavelength (nm) and irradiance
(:math:`W/m^2/nm`)
respectively.  EOM uses set bins and will integrate
those specified in the input file to fit with what is used
in the code.
