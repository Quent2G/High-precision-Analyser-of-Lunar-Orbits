% function [a] = accelHarmonic(r,E,n_max,m_max,Cnm,Snm,GM,R_ref)
% 
% DESCRIPTION
% --------------------------------------------------------------------------
% ACCELHARMONIC provides the inertial acceleration in cartesian coordinates
% due to the gravity field of the central planet. The asymmetry of the 
% gravity field is taken into account by the max degree and max order of
% the harmonics coefficients.
%
% INPUT
% --------------------------------------------------------------------------
% r      = position of the satellite in an inertial frame. 
%          DIMENSION = [km],            SIZE = [3x1]
% E      = rotation matrix from inertial frame to the body fixed frame.
%                                       SIZE = [3x3]
% n      = max degree of the harmonics coefficients.
%                                       SIZE = [1x1]
% m      = max  order of the harmonics coefficients. Set to 0 for taking
%          into account only zonal harmonics.
%                                       SIZE = [1x1]
% Cnm    = normalized harmonics coefficients of cos functions. 
%          use "normalizedharmonics.m" to get Cnm coefficients. 
%                                       SIZE = [nxn]
% Snm    = normalized harmonics coefficients of sin functions.
%          use "normalizedharmonics.m" to get Snm coefficients.
%                                       SIZE = [nxn]
% GM     = gravitational parameter of the central planet
%          DIMENSION = [km^3/s^2],      SIZE = [1x1]
% R_ref  = equatorial radius of the central planet
%          DIMENSION = [km],            SIZE = [1x1]
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
% See also prophpop normalizedharmonics
%
function [a] = accelharmonic(r,E,n,m,Cnm,Snm,GM,R_ref)

% Body-fixed position 
r_bf = E*r;

% Auxiliary quantities
d = norm(r_bf);                     % distance
latgc = asin(r_bf(3)/d);
lon = atan2(r_bf(2)/(d*cos(latgc)),r_bf(1)/(d*cos(latgc)));

[pnm, dpnm] = Legendre(n,m,latgc);

dUdr = 0;
dUdlatgc = 0;
dUdlon = 0;
q3 = 0; q2 = q3; q1 = q2;
for n=0:n
    b1 = (-GM/d^2)*(R_ref/d)^n*(n+1);
    b2 =  (GM/d)*(R_ref/d)^n;
    b3 =  (GM/d)*(R_ref/d)^n;
    for m=0:m
        q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
        q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
        q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
    end
    dUdr     = dUdr     + q1*b1;
    dUdlatgc = dUdlatgc + q2*b2;
    dUdlon   = dUdlon   + q3*b3;
    q3 = 0; q2 = q3; q1 = q2;
end

% Body-fixed acceleration
r2xy = r_bf(1)^2+r_bf(2)^2;

ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

a_bf = [ax ay az]';

% Inertial acceleration 
a = E'*a_bf;
end

% fi [rad]
function [pnm, dpnm] = Legendre(n,m,fi)

    pnm = zeros(n+1,m+1);
    dpnm = zeros(n+1,m+1);

    pnm(1,1)=1;
    dpnm(1,1)=0;
    pnm(2,2)=sqrt(3)*cos(fi);
    dpnm(2,2)=-sqrt(3)*sin(fi);
    % diagonal coefficients
    for i=2:n    
        pnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*cos(fi)*pnm(i,i);
    end
    for i=2:n
        dpnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*((cos(fi)*dpnm(i,i))- ...
                      (sin(fi)*pnm(i,i)));
    end
    % horizontal first step coefficients
    for i=1:n
        pnm(i+1,i)= sqrt(2*i+1)*sin(fi)*pnm(i,i);
    end
    for i=1:n
        dpnm(i+1,i)= sqrt(2*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
    end
    % horizontal second step coefficients
    j=0;
    k=2;
    while(1)
        for i=k:n        
            pnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*pnm(i,j+1))...
                -(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-1,j+1)));
        end
        j = j+1;
        k = k+1;
        if (j>m)
            break
        end
    end
    j = 0;
    k = 2;
    while(1)
        for i=k:n        
            dpnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*dpnm(i,j+1))...
                 +(sqrt(2*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm(i-1,j+1)));
        end
        j = j+1;
        k = k+1;
        if (j>m)
            break
        end
    end
end



