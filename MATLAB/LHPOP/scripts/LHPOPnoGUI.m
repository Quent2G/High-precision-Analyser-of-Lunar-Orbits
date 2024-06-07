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

%% Initialisation Variables
    % ==================================================================================================================================================
    ts = cputime;
    Precise = 0;% 0=Reference else=user's choice --> Computing power
    
    if Precise==0
        HarmD = 70; % maximum degree of the harmonics
        HarmO = 70; % maximum order of the harmonics (set 0 for only zonal harmonics)
        isS = 1;    % Do we considere the Sun as a perturbation
        isE = 1;    % Do we considere the Earth as a perturbation
        RC = 1.3;   % reflection coefficient (0 No SRP - 1 black body - 2 total reflection - 1.3 example)
        AC = 0;     % albedo coefficient     (0 No albedo - 0.3 Earth Albedo coefficient)
        isGR = 1;   % Do we use general relativity
    else
        HarmD = 330;% maximum degree of the harmonics
        HarmO = 330;% maximum order of the harmonics (set 0 for only zonal harmonics)
        isS = 1;    % Do we considere the Sun as a perturbation
        isE = 1;    % Do we considere the Earth as a perturbation
        RC = 1.3;   % reflection coefficient (0 No SRP - 1 black body - 2 total reflection - 1.3 example)
        AC = 0.3;     % albedo coefficient     (0 No albedo - 0.3 Earth Albedo coefficient)
        isGR = 1;   % Do we use general relativity
    end

    Asat = 2;   % satellite area perpendicular to sun direction [m^2]
    msat = 1E3; % satellite mass [kg]
    TStart = '2020 Feb 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    TStop =  '2020 Feb 02 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
    TStep = 30; % Time step [s]
    RefS = 'J2000';   % Reference System (J2000/MOON_ME)
    CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
    X = 8.20151464e+02; % Semi Major Axis
    Y = 1.09123677e+03;    % Eccentricity
    Z = 1.26278945e+03;    % Inclination
    VX = 2.24588615e-01;    % Longitude of Ascending Node
    VY = 1.11772284e+00;    % Argument of Periapsis
    VZ  = -1.13910385e+00;    % Mean Anomaly
    SMA = 1800; % Semi Major Axis
    ECC = 0;    % Eccentricity
    INC = 0;    % Inclination
    LAN = 0;    % Longitude of Ascending Node
    AOP = 0;    % Argument of Periapsis
    MA  = 0;    % Mean Anomaly

%% SPICE LIBRARIES
    % ==================================================================================================================================================
    addpath("mice/lib")
    addpath("mice/src/mice")
    addpath([pwd,'/prop']);

    cspice_kclear; 
    metakernelcheck;
    cspice_furnsh('metakernel.tm');

%% FORCE MODEL
    % ==================================================================================================================================================
    orb.frame.integr = 'J2000';   % integration frame (inertial)
    orb.frame.from   = 'MOON_ME'; % frame where the gravity potential is defined
    orb.frame.to     = orb.frame.integr;
    orb.centralPlanet.stringName = 'Moon';
    orb.pointMasses.stringName   = {'Sun','Earth'};
    % Physical Constants
    orb.const.G           = 6.67428e-20;   % Universal gravitational constant                      [km^3 kg^-1 s^-2]
    % Central Planet
    orb.centralPlanet.RE  = 1737.74;       % Radius of the Moon                                  [km]
    orb.centralPlanet.GM  = 4902.7926024;          % Gravitational parameter of the Moon                   [km^3 s^-2]
    % Point Masses
    orb.pointMasses.M(1)  = 1.9884e30*isS;             % Mass of the Sun                                       [kg]
    orb.pointMasses.M(2)  = 5.97218639e24*isE;        % Mass of the Earth                                      [kg]
    orb.pointMasses.GM    = orb.const.G.*orb.pointMasses.M;        % Gravitational parameter of the Third Bodies           [km^3 s^-2]
    % ==================================================================================================================================================
    orb.pointMasses.numb = length(orb.pointMasses.M);

%% Gravity model
    % ==================================================================================================================================================
    orb.prop.harmonics.degree   = HarmD; % maximum degree of the harmonics
    orb.prop.harmonics.order    = HarmO; % maximum order of the harmonics (set 0 for only zonal harmonics)
    % orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon165x165.txt'];
    orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon_AIUB-GRL350B.txt'];
    % orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon_gggrx_1200a.txt'];
    % ==================================================================================================================================================

    % Harmonics coefficients
    [orb.prop.harmonics.Cnm,orb.prop.harmonics.Snm] = normalizedharmonics(orb.prop.harmonics.filepath,orb.prop.harmonics.degree);

%% Solar Radiation Pressure and Earth Albedo
    % ==================================================================================================================================================
    orb.const.c    = 299792.458; % light speed in the vacuum [km/s];
    orb.const.Ls   = 3.839e26;   % Sun brightness power [W]
    orb.sat.srp.A  = Asat;
    orb.sat.srp.m  = msat;
    orb.sat.srp.CR = RC;
    orb.sat.alb.CR = AC;
    % ==================================================================================================================================================

%% General Relativity
    % ==================================================================================================================================================
    orb.sat.rel    = isGR;
    % ==================================================================================================================================================

%% Analisys start time and stop time
    % ==================================================================================================================================================
    orb.epoch.start = TStart;
    orb.epoch.stop  = TStop;
    % ==================================================================================================================================================
    orb.epoch.str   = char(orb.epoch.start,orb.epoch.stop);   
    orb.epoch.et    = cspice_str2et(orb.epoch.str); % seconds past the Julian Date
    
    orb.epoch.span = orb.epoch.et(1):TStep:orb.epoch.et(end);%+params.orb.periodTrig;
    % ==================================================================================================================================================

%% Initialization
% position and velocity
    % =========================================================================% 
    orb.frame.initstate = RefS;
    if strcmp(CoordT,'Cartesian')
        orb.sat.X0 = [X,Y,Z,VX,VY,VZ];
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_sxform(orb.frame.initstate,orb.frame.to,orb.epoch.et(1));
            orb.sat.X0iner = Rgeog_iner*orb.sat.X0;
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    else
        orb.sat.keplstate(1,1:2) = {'SMA',SMA};        % semi-major axis             [km]
        orb.sat.keplstate(2,1:2) = {'ECC',ECC};        % eccentricity
        orb.sat.keplstate(3,1:2) = {'INC',INC*(pi/180)};  % inclination                 [rad]
        orb.sat.keplstate(4,1:2) = {'LAN',LAN*(pi/180)};   % longitude of ascending node [rad]
        orb.sat.keplstate(5,1:2) = {'AOP',AOP*(pi/180)};   % argument of pericenter      [rad]
        orb.sat.keplstate(6,1:2)  = {'MA',MA*(pi/180)};
        orb.sat.X0 = (cspice_conics([ orb.sat.keplstate{1,2}*(1-orb.sat.keplstate{2,2}),orb.sat.keplstate{2,2},orb.sat.keplstate{3,2},...
                                    orb.sat.keplstate{4,2},orb.sat.keplstate{5,2},orb.sat.keplstate{6,2},...
                                    orb.epoch.et(1),orb.centralPlanet.GM]',orb.epoch.et(1)))';
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_pxform(orb.frame.initstate,orb.frame.to,orb.epoch.et(1));
            orb.sat.X0iner(1:3,1) = Rgeog_iner*orb.sat.X0(1,1:3)';
            orb.sat.X0iner(4:6,1) = Rgeog_iner*orb.sat.X0(1,4:6)';
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    end
    % =========================================================================% 

%%  Propagation of the true state
    options = odeset('RelTol',1e-7,'AbsTol',1e-9);
    [orb.t,orb.XJ2000] = ode45(@prophpop,orb.epoch.span,orb.sat.X0iner,options,orb);
    save('output/ORBdata','orb');
    % ORBIT3D;
    disp(num2str(cputime - ts)+"s")
    
    rmpath([pwd,'/prop']);
    rmpath("mice/lib")
    rmpath("mice/src/mice")