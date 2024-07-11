function f = halo_objective(x)

% halo orbit objective functions

% input

%  x(1) = current normalized x-component of position vector
%  x(2) = current normalized z-component of velocity vector

% output

%  f(1) = normalized x-component of current velocity vector
%  f(2) = normalized z-component of current velocity vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global yi tevent yevent

% set up options for ode45

options = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12, 'Events', @halo_event);

% current position and velocity vectors

ri = yi(1:3);

vi = yi(4:6);

% set normalized x-component of position to current value

vi(1) = x(1);

% set normalized z-component of velocity to current value

vi(3) = x(2);

% maximum normalized search duration

tend = 1000.0;

% solve for y-axis crossing (r_y = 0)

[~, ~, tevent, yevent, ~] = ode45(@halo_eqms, [0.0 tend], [ri vi], options);

% ----------------------------------------------
% current objective functions at y-axis crossing
% ----------------------------------------------

% normalized x-component of current velocity vector

f(1) = yevent(4);

% normalized z-component of current velocity vector

f(2) = yevent(6);



