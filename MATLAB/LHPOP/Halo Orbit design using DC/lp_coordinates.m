function [xl2, yl2] = lp_coordinates

% compute normalized coordinates of the L2 libration point

% output

%  xl2 = normalized x coordinate of L2 libration point
%  yl2 = normalized y coordinate of L2 libration point

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ilp

% L2 libration point (normalized)

ilp = 2.0;

xr1 = -2.0;

xr2 = +2.0;

rtol = 1.0e-12;

[xl2, ~] = brent('clp_func', xr1, xr2, rtol);

yl2 = 0.0;

