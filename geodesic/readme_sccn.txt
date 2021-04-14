This is a modification of geodesic implementation by Danil Kirsanov to
make it compile as a MATLAB MEX file rather than a shared object. The
toolbox .m files are modified to call the mex file.

The original C++ source code is available at
https://code.google.com/archive/p/geodesic/ and the Matlab toolbox can
be downloaded from MathWorks File Exchange.

Compilation Instructions:

  $ mex geodesic.cpp strlcpy.c

requires MATLAB and the boost library to be installed.

Zeynep Akalin Acar, April 2021
