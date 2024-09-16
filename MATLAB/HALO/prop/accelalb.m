% function [accalb] = accelalb(X,sat,const,stringBody,t,stringFrame,stringCentrBody,model)
% 
% DESCRIPTION
% --------------------------------------------------------------------------
% ACCELALB provides the acceleration due to the albedo from a planet on
% a spherical shaped satellite (the cross section is always a circle). 
% The formulation is simple and takes into account only specular reflections.
%
% INPUT 
% --------------------------------------------------------------------------
% X      = current position vector of the satellite [km]
% sat.srp.A  = cross section [m2]
% sat.srp.m  = satellite mass [kg]
% sat.srp.CR = reflection coefficient from the direct radiation (1 black body - 0 no reflection - 2 total reflection)
% sat.alb.CR = reflection coefficient for the albedo (1 black body - 0 no reflection - 2 total reflection)
% const.Ls = Sun brightness power [W]
% const.c  = light speed [km/s]
% stringBody = radiation source (generally SUN)
% t          = seconds past Julian date [s]
% stringFrame = reference frame
% stringCentrBody = central body
% model.centralPlanet.RE = equatorial radius of the central body [km]
%
% OUTPUT 
% --------------------------------------------------------------------------
% accalb = albedo radiation acceleration in the frame specified [km/s^2]
%
% AUTHOR
% --------------------------------------------------------------------------
% Ennio Condoleo,
% Jan 02, 2017 - Rome
% ennio.condoleo@uniroma1.it
%
% See also prophpop accelsrp 
%
function [accalb] = accelalb(X,sat,const,stringBody,t,stringFrame,stringCentrBody,model)

    XB = cspice_spkezr(stringBody,t,stringFrame,'NONE',stringCentrBody);
    XS = cspice_spkezr('SUN',t,stringFrame,'NONE',stringCentrBody);
    
    u = XB(1:3)-X(1:3);
    u = u./norm(u);
    A = sat.srp.A;
    m = sat.srp.m;
    CR = sat.srp.CR;
    c  = const.c*1e3; %m/s
    Ls = const.Ls;
    Calb = sat.alb.CR;
    
    rS = norm(XS(1:3)-XB(1:3));
    r  = norm(XB(1:3));
    rL = norm(X(1:3));
    d  = sqrt(r^2-model.centralPlanet.RE^2);
    dL = norm(X(1:3)-XB(1:3));
    cos_eta  = d/r;
    cos_etaL = ((r^2+dL^2-rL^2)/(2*r*dL));
    if (cos_etaL>cos_eta)&&(dL>d)
        accalb = 0;
    else
        PE = Ls/(4*pi*c*rS^2*1e6);
        PE_mean = PE/4; % solar flux distributed on Earth surface
        P = PE_mean*Calb*(6371/dL)^2;
        accalb = -(P*A*CR/m).*u*1e-3; %km/s^2
    end
end
