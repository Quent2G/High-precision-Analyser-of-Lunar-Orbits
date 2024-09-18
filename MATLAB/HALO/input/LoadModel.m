function [orb] = LoadModel(modeldef,orb)
%LOADMODEL The model used by the propagator is defined here

    if modeldef == "Ref"
        HarmD = 150; % e-3 , maximum degree of the harmonics
        HarmO = 150; % e-3 , maximum order of the harmonics (set 0 for only zonal harmonics)
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

    elseif modeldef == "Fast"
        HarmD = 8; % e-3 , maximum degree of the harmonics
        HarmO = 8; % e-3 , maximum order of the harmonics (set 0 for only zonal harmonics)
        isE = 1;    % e-8 , Do we considere the Earth as a perturbation
        isS = 0;    % e-11, Do we considere the Sun as a perturbation
        RC = 1.3;   % e-11, reflection coefficient (0 No SRP - 1 black body - 2 total reflection - 1.3 example)
        isGR = 0;   % e-13, Do we use general relativity
        isJ = 0;    % e-14, Do we considere Jupiter as a perturbation
        HarmDE = 0; % e-14, maximum degree of the Earth harmonics
        HarmOE = 0; % e-14, maximum order of the Earth harmonics (set 0 for only zonal harmonics)
        AC = 0;   % e-16, albedo coefficient     (0 No albedo - 0.3 Earth Albedo coefficient)
        
        Asat = 2;   % satellite area perpendicular to sun direction [m^2]
        msat = 1E3; % satellite mass [kg]
    end
    
%% FORCE MODEL
    % ==================================================================================================================================================
    orb.frame.integr = 'J2000';   % integration frame (inertial)
    orb.frame.from   = 'MOON_PA'; % frame where the lunar gravity potential is defined
    orb.frame.fromE   = 'ITRF93'; % frame where the Earth gravity potential is defined
    orb.frame.to     = orb.frame.integr;
    orb.centralPlanet.stringName = 'Moon';
    orb.Earth.stringName = 'Earth';
    orb.pointMasses.stringName   = {'Sun','JUPITER BARYCENTER'};
    % Physical Constants
    orb.const.G           = 6.67428e-20;   % Universal gravitational constant                      [km^3 kg^-1 s^-2]
    % Moon
    orb.centralPlanet.RE  = 1738;       % Radius of the Moon                                  [km]
    orb.centralPlanet.GM  = 4902.7999671;  % Gravitational parameter of the Moon                   [km^3 s^-2]
    % Earth
    orb.Earth.RE  = 6378.1363;               % Radius of the Earth                                  [km]
    orb.Earth.GM  = 398600.4415*isE;          % Gravitational parameter of the Earth                   [km^3 s^-2]
    % Point Masses : Sun and Jupiter
    orb.pointMasses.M(1)  = 1.9884e30*isS;             % Mass of the Sun                                       [kg]
    orb.pointMasses.M(2)  = 1.89813e27*isJ;             % Mass of the Sun                                       [kg]
    orb.pointMasses.GM    = orb.const.G.*orb.pointMasses.M;        % Gravitational parameter of the Third Bodies           [km^3 s^-2]
    orb.pointMasses.numb = length(orb.pointMasses.M);
    % ==================================================================================================================================================
    
%% Gravity models
    % ==================================================================================================================================================
    orb.prop.harmonics.degree   = HarmD; % maximum degree of the harmonics
    orb.prop.harmonics.order    = HarmO; % maximum order of the harmonics (set 0 for only zonal harmonics)
    orb.prop.harmonics.filepath = [cd,'/input/gravity_models/Moon_AIUB-GRL350B.txt'];

    orb.prop.harmonics.degreeE   = HarmDE; % maximum degree of the harmonics
    orb.prop.harmonics.orderE    = HarmOE; % maximum order of the harmonics (set 0 for only zonal harmonics)
    orb.prop.harmonics.filepathE = [cd,'/input/gravity_models/Earth_EGM2008.txt'];
   
    % Harmonics coefficients
    [orb.prop.harmonics.Cnm,orb.prop.harmonics.Snm] = normalizedharmonics(orb.prop.harmonics.filepath,orb.prop.harmonics.degree);
    [orb.prop.harmonics.ECnm,orb.prop.harmonics.ESnm] = normalizedharmonics(orb.prop.harmonics.filepathE,orb.prop.harmonics.degreeE);
    % ==================================================================================================================================================

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

end

