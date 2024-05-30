% function [accsrp] = accelsrp(X,srp,const,stringBody,t,stringFrame,stringCentrBody,model)
% 
% DESCRIPTION
% --------------------------------------------------------------------------
% ACCELSRP provides the acceleration due to the solar radiation pressure acting on
% a spherical shaped satellite (the cross section is always a circle). 
% The formulation is simple and takes into account only specular reflections.
%
% INPUT 
% --------------------------------------------------------------------------
% X      = current position vector of the satellite [km]
% srp    = structure with fields
%           A  = cross section [m2]
%           m  = satellite mass [kg]
%           CR = reflection coefficient (1 black body - 0 no reflection - 2 total reflection)
% const  = structure with fields
%           Ls = Sun brightness power [W]
%           c  = light speed [km/s]
% stringBody = radiation source (generally SUN)
% t          = seconds past Julian date [s]
% stringFrame = reference frame
% stringCentrBody = central body
% model.centralPlanet.RE = equatorial radius of the central body [km]
%
% OUTPUT 
% --------------------------------------------------------------------------
% accsrp = solar radiation pressure acceleration in the frame specified [km/s^2]
%
% AUTHOR
% --------------------------------------------------------------------------
% Ennio Condoleo,
% Jan 02, 2017 - Rome
% ennio.condoleo@uniroma1.it
%
% See also prophpop accelalb
%

function [accsrp] = accelsrp(X,srp,const,stringBody,t,stringFrame,stringCentrBody,model)

    XB = cspice_spkezr(stringBody,t,stringFrame,'NONE',stringCentrBody);
    u = XB(1:3)-X(1:3);
    u = u./norm(u);
    A = srp.A;
    m = srp.m;
    CR = srp.CR;
    c  = const.c*1e3; %m/s
    Ls = const.Ls;
    
    r  = norm(XB(1:3));
    rL = norm(X(1:3));
    d  = sqrt(r^2-model.centralPlanet.RE^2);
    dL = norm(X(1:3)-XB(1:3));
    cos_eta  = d/r;
    cos_etaL = ((r^2+dL^2-rL^2)/(2*r*dL));
    if (cos_etaL>cos_eta)&&(dL>d)
        accsrp = 0;
    else
        P = Ls/(4*pi*c*dL^2*1e6);
        accsrp = -(P*A*CR*1e-3/m).*u; %km/s^2
    end
end
