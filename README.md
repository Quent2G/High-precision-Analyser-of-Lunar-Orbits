## Low Lunar Orbit Propagator

## Description

The purpose of this project is to improve an existing lunar orbit propagator in order to make it accurate on low lunar orbit, and later on DRO and NRHO. Thus the propagator will allow us to find frozen orbits and to properly simulate the evolution of a spacecraft in the vicinity of the moon.

## Initialisation

After copying the project on your local machine (git clone [url]), you have to:

1. (If you are a Windows user, skip this part (the windows mice is already installed), or else:) Download the mice folder depending on your machine OS at this link : [mice](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html)  
   Relace the mice folder in the MATLAB/LHPOP folder with your one.
2. Download the de441.bsp at [JPL de441](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/) --> no changes needed  
    or the de430.bsp ephemeris at [NAIF de430](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/) --> need to change the metakernelcheck.m  
   Place it in MATLAB/LHPOP/ker.

Attention for Mac users with a Mac processor (no problem with Intel processors), you may need to recompile the mexfile following these instructions for Mac Silicon processors, e.g., M2:

1. Download the files mice.tar.Z and importMice.csh from NAIF(https://naif.jpl.nasa.gov/naif/toolkit_MATLAB_MacIntel_OSX_AppleC_MATLAB9.x_64bit.html) into a directory, such as Build;
2. Open the Build directory in Matlab, and in the Matlab Command Window, run:
   a. system('/bin/csh -f importMice.csh'), followed by,  
   b. system('/bin/csh -f makeall.csh');
3. A new mex file mice.mexmaci64 should be generated in the lib sub directory;
4. Copy the whole mice directory to your project directory.

## Usage

The MATLAB folder contains all you need to run the Matlab propagator.  
The Python folder contains all you need for the data visualisation python files.  
The Docs folder contains every other document.

## Authors and acknowledgment

Creator of the Github and first contributor to the project: Quentin Granier, 3 months Research Associate.  
Supervisor of the Project : Dr Yang Yang.
