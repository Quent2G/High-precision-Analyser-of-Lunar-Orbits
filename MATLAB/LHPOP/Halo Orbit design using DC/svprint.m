function svprint(r, v)

% print normalized position and velocity vectors and magnitudes

% input

%  r = normalized position vector
%  v = normalized velocity vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized position magnitude

rmag = norm(r);

% normalized velocity magnitude

vmag = norm(v);

% print normalized state vector and magnitudes

fprintf ('\n          rx                     ry                     rz                    rmag');

fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e  \n', r(1), r(2), r(3), rmag);

fprintf ('\n          vx                     vy                     vz                    vmag');

fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e  \n\n', v(1), v(2), v(3), vmag);


