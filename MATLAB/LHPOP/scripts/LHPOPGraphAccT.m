%% Initialisation Variables
    % ==================================================================================================================================================
    ts = cputime;
    
    HarmD = 70; % e-3 , maximum degree of the harmonics
    HarmO = 70; % e-3 , maximum order of the harmonics (set 0 for only zonal harmonics)
    isE = 1;    % e-8 , Do we considere the Earth as a perturbation
    isS = 1;    % e-11, Do we considere the Sun as a perturbation
    RC = 1.3;   % e-11, reflection coefficient (0 No SRP - 1 black body - 2 total reflection - 1.3 example)
    isGR = 1;   % e-13, Do we use general relativity
    isJ = 1;    % e-14, Do we considere Jupiter as a perturbation
    HarmDE = 3; % e-14, maximum degree of the Earth harmonics
    HarmOE = 3; % e-14, maximum order of the Earth harmonics (set 0 for only zonal harmonics)
    AC = 0.3;   % e-16, albedo coefficient     (0 No albedo - 0.3 Earth Albedo coefficient)

    Asat = 2;   % satellite area perpendicular to sun direction [m^2]
    msat = 1E3; % satellite mass [kg]
    TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    TStop =  '2022 Dec 10 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
    TStep = 30; % Time step [s]
    RefS = 'J2000';   % Reference System (J2000/MOON_ME)
    CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
    X = -1642.1842792591378; % Semi Major Axis
    Y = 62149.566213382524;    % Eccentricity
    Z = 32159.59842593645;    % Inclination
    VX = 0.3549830252637374;    % Longitude of Ascending Node
    VY = 0.017074626095640213;    % Argument of Periapsis
    VZ  = -0.016135387025700578;    % Mean Anomaly

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
    orb.frame.fromE   = 'IAU_EARTH'; % frame where the gravity potential is defined
    orb.frame.to     = orb.frame.integr;
    orb.centralPlanet.stringName = 'Moon';
    orb.Earth.stringName = 'Earth';
    orb.pointMasses.stringName   = {'Sun','JUPITER BARYCENTER'};
    % Physical Constants
    orb.const.G           = 6.67428e-20;   % Universal gravitational constant                      [km^3 kg^-1 s^-2]
    % Moon
    orb.centralPlanet.RE  = 1737.74;       % Radius of the Moon                                  [km]
    orb.centralPlanet.GM  = 4902.7926024;  % Gravitational parameter of the Moon                   [km^3 s^-2]
    % Earth
    orb.Earth.RE  = 6371.00;               % Radius of the Earth                                  [km]
    orb.Earth.M  = 5.97218639e24*isE;        % Mass of the Earth                                      [kg]
    orb.Earth.GM  = orb.const.G.*orb.Earth.M;          % Gravitational parameter of the Earth                   [km^3 s^-2]
    % Point Masses : Sun and Jupiter
    orb.pointMasses.M(1)  = 1.9884e30*isS;             % Mass of the Sun                                       [kg]
    orb.pointMasses.M(2)  = 1.89813e27*isJ;             % Mass of the Sun                                       [kg]
    orb.pointMasses.GM    = orb.const.G.*orb.pointMasses.M;        % Gravitational parameter of the Third Bodies           [km^3 s^-2]
    % ==================================================================================================================================================
    orb.pointMasses.numb = length(orb.pointMasses.M);

%% Gravity models
    % ==================================================================================================================================================
    orb.prop.harmonics.degree   = HarmD; % maximum degree of the harmonics
    orb.prop.harmonics.order    = HarmO; % maximum order of the harmonics (set 0 for only zonal harmonics)
    % orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon165x165.txt'];
    orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon_AIUB-GRL350B.txt'];
    % orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon_gggrx_1200a.txt'];

    orb.prop.harmonics.degreeE   = HarmDE; % maximum degree of the harmonics
    orb.prop.harmonics.orderE    = HarmOE; % maximum order of the harmonics (set 0 for only zonal harmonics)
    orb.prop.harmonics.filepathE = [cd,'/input/gravity_models/Earth_C-EGM2008.txt'];
    % ==================================================================================================================================================

    % Harmonics coefficients
    [orb.prop.harmonics.Cnm,orb.prop.harmonics.Snm] = normalizedharmonics(orb.prop.harmonics.filepath,orb.prop.harmonics.degree);
    [orb.prop.harmonics.ECnm,orb.prop.harmonics.ESnm] = normalizedharmonics(orb.prop.harmonics.filepathE,orb.prop.harmonics.degreeE);

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
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    
    global TextAcc;
    TextAcc = fopen('output/acc.txt',"w");
    [orb.t,orb.XJ2000] = ode45(@prophpopGraphAccT,orb.epoch.span,orb.sat.X0iner,options,orb);
    fclose(TextAcc);
    save('output/ORBdata','orb');
    % ORBIT3D;
    disp(num2str(cputime - ts)+"s")
    
    rmpath([pwd,'/prop']);
    rmpath("mice/lib")
    rmpath("mice/src/mice")