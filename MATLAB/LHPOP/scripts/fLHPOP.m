% LHPOP is a High Precision Orbit Propagator here built for the propagation of a lunar orbiter.
% The motion is integrated under 
% - the Keplerian attraction force
% - the perturbations due to the asymmetry of the lunar gravity field (up to 165x165 harmonics)
% - the third body attractions of the Sun and Earth
% - the solar radiation pressure
% - the Earth's albedo
% - the General relativity
%
% When you run the code, a GUI appears and the following inputs are required:
% - degree of the lunar gravity field
% - order  of the lunar gravity field
% - Sun's   perturbation
% - Earth's perturbation
% - General relativity
% - Solar radiation pressure
% - mass of the satellite 
% - cross section of the satellite (spherical volume)
% - reflection coefficient for the solar radiation
% - reflection coefficient for the Earth's albedo
% - initial state of the orbiter (keplerian/cartesian)
% - reference frame
% - start time
% - stop time
% - step/integration time
%     
% The integration is performed in the "J2000" inertial frame, and the output is the state and time vectors.
% AT the end of the integration the code "ORBIT3D" is launched, and after the loading of the orbiter
% data, the trajectory is plotted.
%
% IMPORTANT: it is required to install mice toolbox from
% https://naif.jpl.nasa.gov/naif/toolkit_MATLAB_PC_Windows_VisualC_MATLAB8.x_64bit.html
%
% Note: work is in progress to improve the code. 
%
% (c)Copyright 2017 - Ennio Condoleo
% Sapienza University of Rome, 08-25-2017
% ennio.condoleo@uniroma1.it

    ts = cputime;

%% SPICE LIBRARIES
    addpath("mice/lib");
    addpath("mice/src/mice");
    addpath('prop');
    addpath('input');
    
    cspice_kclear; 
    metakernelcheck;
    cspice_furnsh('metakernel.tm');

%% FORCE MODEL
    orb = struct();
    orb = LoadModel("Ref",orb);
    
%% Initialization
    orb = LoadState("DRO14l",orb);

%%  Propagation of the true state
    
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    [orb.t,orb.XJ2000] = ode45(@prophpop,orb.epoch.span,orb.sat.X0iner,options,orb);

    % X = fsolveConv(orb.sat.X0iner,orb);
    % [~,orb.XC] = ode45(@prophpop,orb.epoch.span,X,options,orb);

    save('output/ORBdata','orb');
    % ORBIT3D;
    disp(num2str(cputime - ts)+"s")
    
%% Close Path
    rmpath('prop');
    rmpath('input');
    rmpath("mice/lib")
    rmpath("mice/src/mice")

