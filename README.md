# High-precision Analyser of Lunar Orbits

# Description

HALO is a mission design tool build on a high accuracy propagator for lunar orbits. It intends to help to study any orbits in the vicinity of the Moon. HALO comes with a scientific paper that precisely define all the models and algorithms used and go through a sensitivity analysis on lunar orbits of interest. This paper High-precision_Analyser_of_Lunar_Orbit.pdf can be found at the root of the GitHub, it is highly recommended to get a good look at the paper when using HALO.

# Initialisation

After copying the project on your local machine (git clone [url]), you have to:

1. If you are a Windows user, skip this part (the windows mice is already installed), or else: Download the mice folder depending on your machine OS at this link : [mice](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html)  
   Replace the mice folder in the MATLAB/HALO folder with your one.
2. Download the de430.bsp ephemeris at [NAIF de430](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/).   
   Place it in MATLAB/HALO/ker and in Python/input.

Attention for Mac users with a Mac processor (no problem with Intel processors), you may need to recompile the mexfile following these instructions for Mac Silicon processors, e.g., M2:

1. Download the files mice.tar.Z and importMice.csh from [NAIF](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB_MacIntel_OSX_AppleC_MATLAB9.x_64bit.html) into a directory, such as Build;
2. Open the Build directory in Matlab, and in the Matlab Command Window, run:
   - system('/bin/csh -f importMice.csh');
   - system('/bin/csh -f makeall.csh');
3. A new mex file mice.mexmaci64 should be generated in the lib sub directory;
4. Copy the whole mice directory to your project directory in MATLAB/HALO.

### List of all Python modules
- matplotlib
- spiceypy
- numpy
- tqdm
- pandas
- pymatreader

# Usage

The MATLAB folder contains all you need to run the tool.  
The Python folder contains all you need for the data visualisation.

The main usage of HALO is to construct a mission by changing the input files like LoadSequential.m on MATLAB, then run the mission by running the main script HALO.m, and finally use python scripts like Visualisation.py to process the output file. Because of relative path, HALO.m as to be run from the HALO folder. This was chosen to avoid windows path length limitation.   
Read the following parts to learn what can be done with HALO.

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

### Process an ephemeris 3-body-problem orbit
You can get a CR3BP orbit from the [JPL website](https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro). This orbit is a solution in a simplified model and in the Earth-Moon rotational frame.   
CR3BPfitting.py will fit the simplified model on an ephemeris one at a specific epoch (chosen at 1 Jan 2024 in the example) and convert the initial position to be used for propagation.   
Then a sequence of type "TBPOtim" allows to compute an ephemeris orbit from this fitted initial guess.    
The following sequential gives an example for a NRHO close to the Capstone one:  

   orb = LoadState("NRHO",orb);  
   orb.seq.Time = cspice_str2et(orb.sat.t0);    
   orb.seq.a.type = "TBPOptim";  
   orb.seq.a.T = 5.7444e+05;  

Finally, Visualisation.py allows to visualise the initial guess and the converged solution by putting the Converged option to 1. The RotationalF option can also be put to 1 to observe the closed orbit in the rotational frame.   

The same can be done for DRO14, with a period of 14 days, but the convergence is a bit harder so that we have to change n=2 (in the algorithm's file) and use only two periods due to the long period of the DRO.

### Some useful functions
sp.dafec(sp.spklef("input/LRO_ES_90_202003_GRGM900C_L600.bsp"),1000)    :    Allows to assess the comment of a bsp file
sp.datetime2et(sp.datetime(y, m, d, h, m, s))                           :    Converts a UTC time to ET

See also part 4.3.5 of the paper for time and coordinate systems handling in either MATLAB and Python.

# Authors and acknowledgment

Creator of the Github and first contributor to the project: Quentin Granier, Research Associate.  
Supervisor of the Project : Dr Yang Yang, contact: yiyinfeixiong@gmail.com.
