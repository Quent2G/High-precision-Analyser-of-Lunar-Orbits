% function [a] = accelpntmasses(r,pointMasses,GMpointMasses,t,frame,centralPlanet)
%
% DESCRIPTION
% --------------------------------------------------------------------------
% ACCELPNTMASSES provides the inertial acceleration in cartesian coordinates
% due to the gravity field of third bodies. The ephemerides are given by
% spice libraries (DE430).
%
% INPUT
% --------------------------------------------------------------------------
% r             = position of the satellite in an inertial frame. 
%                   DIMENSION = [km],           SIZE = [3x1]
% pointMasses   = cell with name string of the third body/ies.
%                                               SIZE = [1xN]
% GMpointMasses = gravitation parametr of the third body/ies.
%                   DIMENSION = [km^3/s^2],     SIZE = [1xN]
% t             = seconds past 2000 Jan 01 12:00:00.000
%                   DIMENSION = [s],     SIZE = [1x1]
% frame         = inertial frame string (e.g. 'J2000','ICRF') 
% centralPlanet = central planet string;
%
% OUTPUT
% --------------------------------------------------------------------------
% a      = acceleration of the satellite in an inertial frame. 
%          DIMENSION = [km/s^2],            SIZE = [3x1]
%
% AUTHOR
% --------------------------------------------------------------------------
% Ennio Condoleo,
% Jan 02, 2017 - Rome
% ennio.condoleo@uniroma1.it
%
% See also prophpop 
%
function [acc_pointMasses] = accelpntmasses(r,pointMasses,GM,t,frame,centralPlanet,model) %#ok<INUSL>

    acc_pointMasses = zeros(3,1);
    for j = 1:size(pointMasses,2)
        temp = cspice_spkezr(pointMasses{j},t,frame,'NONE',centralPlanet); %#ok<NASGU>
        jstr = num2str(j);       
        eval(['r0',jstr,' = temp(1:3);']);
        eval(['r',jstr,'S = r-r0',jstr,';']);
        eval(['r0',jstr,'_3 = norm(r0',jstr,')^3;']);
        eval(['r',jstr,'S_3 = norm(r',jstr,'S)^3;']);
        eval(['acc_pointMasses = acc_pointMasses + GM(j).*(-r0',jstr,...
            './r0',jstr,'_3 - r',jstr,'S./r',jstr,'S_3);']);
    end
end
