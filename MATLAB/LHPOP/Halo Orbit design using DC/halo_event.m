function [value, isterminal, direction] = halo_event(~, y)

% x-axis crossing event function

% input

%  y = normalized state vector

% output

%  value = normalized y-component of position

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized y position component

value = y(2);

isterminal = 1;

direction =  -1;



