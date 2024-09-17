# High-precision Analyser of Lunar Orbits

# Description

HALO is a mission design tool build on a high accuracy propagator for lunar orbits. It intends to help to study any orbits in the vicinity of the Moon. HALO comes with a scientific paper that precisely define all the models and algorithms used and go through a sensitivity analysis on lunar orbits of interest.

# Initialisation

After copying the project on your local machine (git clone [url]), you have to:

1. If you are a Windows user, skip this part (the windows mice is already installed), or else: Download the mice folder depending on your machine OS at this link : [mice](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html)  
   Replace the mice folder in the MATLAB/LHPOP folder with your one.
2. Download the de430.bsp ephemeris at [NAIF de430](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/)
   Place it in MATLAB/LHPOP/ker.

Attention for Mac users with a Mac processor (no problem with Intel processors), you may need to recompile the mexfile following these instructions for Mac Silicon processors, e.g., M2:

1. Download the files mice.tar.Z and importMice.csh from [NAIF](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB_MacIntel_OSX_AppleC_MATLAB9.x_64bit.html) into a directory, such as Build;
2. Open the Build directory in Matlab, and in the Matlab Command Window, run:
   - a. system('/bin/csh -f importMice.csh'), followed by,
   - b. system('/bin/csh -f makeall.csh');
3. A new mex file mice.mexmaci64 should be generated in the lib sub directory;
4. Copy the whole mice directory to your project directory in MATLAB/LHPOP.

# Usage

The MATLAB folder contains all you need to run the tool.  
The Python folder contains all you need for the data visualisation.

The main usage of HALO is to construct a mission by changing the input files like LoadSequential.m on MATLAB, then run the mission by running the main script HALO.m, and finally use python scripts like Visualisation.py to process the output file.

### Architechture of the MATLAB
The 2 metakernel files allow to dialogue with SPICE by loading the relevant kernel stored in the "ker" folder   
HALO.m is the main file that allows to run the tool   
The "prop" folder gathers all the calculation related files with notably prophpop.m defining the motion equation thanks to the other files   
The "mice" folder gathers the SPICE files   
The output folder contains all the information on the processed mission   
The input folder contains among others the mission design informations as the sequential, the model or the initial state.

### Architechture of the Python
Visualisation.py allows to plot a propagation after processing, the user can choose the frame, or choose to add the converged trajectory if there is one, same with a CR3BP trajectory.   
TrajectoryErrors.py allows to compare the propagation of a true ephemeris state of a satellite with its real evolution.   
CR3BPfitting.py converts the initial state from a CR3BP orbit chosen on https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro in the Earth-Moon rotationnal and barycentered frame into the ephemeris Moon centered inertial frame in order to use it as an initial state for the propagator.   
process.py gathers all the processing functions used in the other python files.

# Authors and acknowledgment

Creator of the Github and first contributor to the project: Quentin Granier, Research Associate.  
Supervisor of the Project : Dr Yang Yang, contact: yiyinfeixiong@gmail.com.
