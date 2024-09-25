% function [X] = prophpop(t,X0,model)
% 
% DESCRIPTION
% --------------------------------------------------------------------------
% PROPHPOP integrates the motion of a spacecraft under the
% attraction of third bodies and considering the asymettry of
% the central planet gravity field. It has to be used in the call
% of ode functions as follows
% [t,X] = ode113(@prophpop,tspan,X0,options,model)
%
% INPUT ODE113 FOR PROPHPOP
% --------------------------------------------------------------------------
% tspan  = time vector for the integration [s]
% X0     = initial state [km for position, km/s for velocity] (1x6)
% model  = struct with the following fields
%           model.pointMasses.stringName    = cell array of the name strings of the third bodies
%           model.pointMasses.GM            =      array of the gravitational parameters of the third bodies [km^3/s^2]
%           model.frame.integr              = string of the inertial reference system for the integration ('J2000' suggested)
%           model.frame.from                = string of the body-fixed reference system (e.g. 'MOON_ME','IAU_EARTH','IAU_MARS', ITRF93')
%           model.frame.to                  = string of the inertial reference system for the integration (as  model.frame.integr)
%           model.prop.harmonics.degree     = max degree for the harmonics coefficients (with added E for Earth)
%           model.prop.harmonics.order      = max order  for the harmonics coefficients (with added E for Earth)
%           model.prop.harmonics.Cnm        = cosine functions coefficients up to max degree (use normalizedharmonics.m to get them) (with added E for Earth)
%           model.prop.harmonics.Snm        =   sine functions coefficients up to max degree (use normalizedharmonics.m to get them) (with added E for Earth)
%           model.centralPlanet.stringName  = name string of the central planet
%           model.centralPlanet.GM          = gravitational parameter of the central planet [km^3/s^2]
%           model.centralPlanet.RE          = equatorial radius       of the central planet [km]
%           model.Earth.stringName          = name string of the central planet
%           model.Earth.GM                  = gravitational parameter of the Earth [km^3/s^2]
%           model.Earth.RE                  = equatorial radius       of the Earth [km]
%           model.const.c                   = light speed in the vacuum [km/s]
%           model.const.Ls                  = Sun brightness power [W]
%           model.sat.rel                   = activator for the general relativity perturbation [1/0]
%           model.sat.srp.A                 = cross section [m2]
%           model.sat.srp.m                 = satellite mass [kg]
%           model.sat.srp.CR                = reflection coefficient from the direct radiation (1 black body - 0 no reflection - 2 total reflection)
%           model.sat.alb.CR                = reflection coefficient for the albedo (1 black body - 0 no reflection - 2 total reflection)
%
% OUTPUT FROM ODE
% --------------------------------------------------------------------------
% t      = integration time [s]
% X      = state vector [km,km/s]
%
% AUTHORS
% --------------------------------------------------------------------------
% Ennio Condoleo,
% Jan 02, 2017 - Rome
% ennio.condoleo@uniroma1.it
%
% 2024 - Quentin Granier
% Research project at University of New South Wales
% Supervisor's mail: yiyinfeixiong@gmail.com
%
% See also accelpntmasses accelharmonic accelsrp accelalb
%
function [X] = prophpop(t,X0,model)
    
    % state in J2000 reference system centred at the central planet
    xJ2000 = X0(1:3);
    vJ2000 = X0(4:6);
    
    % acceleration due to the attraction of the Sun
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_pointMasses = accelpntmasses(xJ2000,model.pointMasses.stringName,model.pointMasses.GM,...    
        t,model.frame.integr,model.centralPlanet.stringName);
    % ---------------------------------------------------------------------------------------------------------- %

    % acceleration due to the central planet gravity field
    % ---------------------------------------------------------------------------------------------------------- %
    Rgeog_iner = cspice_pxfrm2(model.frame.from,model.frame.to,t,t);
    acc_centralPlanet = accelharmonic(xJ2000,Rgeog_iner',model.prop.harmonics.degree,model.prop.harmonics.order,...
        model.prop.harmonics.Cnm,model.prop.harmonics.Snm,model.centralPlanet.GM,model.centralPlanet.RE);
    % ---------------------------------------------------------------------------------------------------------- %   

    % acceleration due to the Earth gravity field
    % ---------------------------------------------------------------------------------------------------------- %
    REarth_iner = cspice_pxfrm2(model.frame.fromE,model.frame.to,t,t);
    R_Moon_Earth = cspice_spkezr('EARTH',t,model.frame.to,'NONE','MOON');
    R_Moon_Earth = R_Moon_Earth(1:3);
    acc_earth = accelharmonic(xJ2000-R_Moon_Earth,REarth_iner',model.prop.harmonics.degreeE,model.prop.harmonics.orderE,...
                model.prop.harmonics.ECnm,model.prop.harmonics.ESnm,model.Earth.GM,model.Earth.RE) + ...
                -accelharmonic(-R_Moon_Earth,REarth_iner',model.prop.harmonics.degreeE,model.prop.harmonics.orderE,...
                model.prop.harmonics.ECnm,model.prop.harmonics.ESnm,model.Earth.GM,model.Earth.RE);
    % ---------------------------------------------------------------------------------------------------------- %   
    
    % acceleration due to the general relativity
    % ---------------------------------------------------------------------------------------------------------- %    
    gamma = 1; beta = 1; c = model.const.c;
    acc_genRel = model.sat.rel.*model.centralPlanet.GM/(c^2*norm(xJ2000)^3)*((2*(beta+gamma)*model.centralPlanet.GM/norm(xJ2000)-...
        gamma*norm(vJ2000)^2).*xJ2000+2*(1+gamma)*xJ2000'*vJ2000.*vJ2000);
    % ---------------------------------------------------------------------------------------------------------- %

    % acceleration due to the solar radiation pressure of Sun
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_SRPSun = accelsrp(xJ2000,model.sat.srp,model.const,'SUN',t,model.frame.integr,model.centralPlanet.stringName,model);
    % ---------------------------------------------------------------------------------------------------------- %
    
    % acceleration due to the Earth albedo
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_alb = accelalb(xJ2000,model.sat,model.const,'EARTH',t,model.frame.integr,model.centralPlanet.stringName,model);
    % ---------------------------------------------------------------------------------------------------------- %

    % equations of motion
    % ---------------------------------------------------------------------------------------------------------- %
    xdot   = vJ2000;
    vdot   = acc_centralPlanet + acc_earth + acc_pointMasses + acc_SRPSun + acc_alb + acc_genRel;
    % ---------------------------------------------------------------------------------------------------------- %

    % output
    % ---------------------------------------------------------------------------------------------------------- %
    X = [xdot;vdot];
    % ---------------------------------------------------------------------------------------------------------- %   
end
