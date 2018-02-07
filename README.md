# CosmoLSSwrapper
Python wrapper for the CosmoLSS likelihood. It can be used as a
template for wrapping other CosmoMC likelihoods. 

# Installation
- Download the install script ````CosmoLSSinstaller.sh```` to you
MontePython root folder.
- If your C and Fortran compiler are different from gcc, change
CC=gcc and FC=gcc in the last line of the install script. Example:
CC=gcc-6 FC=gcc-6 make
- Run ````bash  CosmoLSSinstaller.sh````. Beware that it will overwrite
the CosmoLSS directory in montepython/likelihoods if it exist and also
the CosmoLSS directory in data.

# Some details
The CosmoLSS likelihood depends mainly on three things:
1) Several helper classes defined in CosmoMC
2) a set of functions provided by CAMB, like cosmological
distances.
3) data that are computed by CAMB but made available through the
CosmoMC class ````TCosmoTheoryPredictions```` like the matter power
spectrum.

Item 1 is handled by using part of Antony Lewis [forutils] (https://github.com/cmbant/forutils)
package, which is a subset of CosmoMC for general purpose.
Item 2 and 3 are taken care of by creating Fortran interpolation
objects of all the relevant data from CLASS.

The cosmological parameters in CosmoMC are stored in the class
````CMBParams````, which in this implementation became a derived type
which must be updated every time the Python likelihood is called.

In the future, if all likelihoods restricted themselves to using only
the two classes ````CMBParams```` and
````TCosmoTheoryPredictions````, it might be worth wrapping these two
classes in a more complete way.
