function dx = halo_eqms(~, y)

% first order circular restricted three-body equations of motion

% input

%  y(1) = x-component of normalized position
%  y(2) = y-component of normalized position
%  y(3) = z-component of normalized position
%  y(4) = x-component of normalized velocity
%  y(5) = y-component of normalized velocity
%  y(6) = z-component of normalized velocity

% output

%  dx(1:6) = normalized acceleration vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu

r1 = zeros(3, 1);

r2 = zeros(3, 1);

% extract normalized position vector

r = y(1:3);

% extract normalized velocity vector

v = y(4:6);

% normalized location of minor body

rb1(1) = -mu;
rb1(2) = 0.0;
rb1(3) = 0.0;

% normalized location of major body

rb2(1) = 1.0 - mu;
rb2(2) = 0.0;
rb2(3) = 0.0;

for i = 1:1:3
    
    % normalized position vector from minor body to spacecraft
    
    r1(i) = r(i) - rb1(i);
    
    % normalized position vector from major body to spacecraft
    
    r2(i) = r(i) - rb2(i);
    
end

r13 = norm(r1)^3;

r23 = norm(r2)^3;

c1  = (1.0 - mu) / r13;

c2  = mu / r23;

% -----------------------------
% normalized integration vector
% -----------------------------

dx = [v(1)
    
v(2)

v(3)

2.0 * v(2) + r(1) - c1 * (r(1) + mu) - c2 * (r(1) - 1.0 + mu)

-2.0 * v(1) + r(2) - c1 * r(2) - c2 * r(2)

-c1 * r(3) - c2 * r(3)];

end
