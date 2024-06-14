## Low Lunar Orbit Propagator

## Description
The purpose of this project is to improve an existing lunar orbit propagator in order to make it accurate on low lunar orbit, and later on DRO and NRHO. Thus the propagator will allow us to find frozen orbits and to properly simulate the evolution of a spacecraft in the vicinity of the moon.

## Initialisation
After copying the project on your local machine (git clone [url]), you have to:  
  1) Download the mice folder depending on your machine OS at this link : [mice](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html)  
Place the mice folder directly in the LHPOP folder  
  2) Download the de430.bsp ephemeris at https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/  
Place it in MATLAB/LHPOP/ker  
  
Attention for Mac users with a mac processor (no problem with Intel processors), you may need to recompile the mexfile following these instructions:  
--> write instructions here <--

## Usage
The MATLAB folder contains all there need to run the Matlab propagator.  
The Python folder contains all there need for the data visualisation python files.  
The Docs folder contains every other document.

## Authors and acknowledgment
Creator of the Github and first contributor to the project: Quentin Granier, 3 months Research Associate  
Supervisor of the Project : Dr Yang Yang
