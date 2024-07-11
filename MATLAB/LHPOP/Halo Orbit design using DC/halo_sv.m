function [r_crtbp, v_crtbp, period] = halo_sv(amp_z, t)

% halo orbit normalized state vector and period

% Richardson's analytic solution

% input

%  amp_z = halo z amplitude (kilometers)
%  t     = current normalized time

% output

%  r      = normalized halo orbit position vector
%  period = normalized halo orbit period

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu hclass dist ilp

r_crtbp = zeros(3, 1);

v_crtbp = zeros(3, 1);

ilp = 2;

% compute L2 libration point normalized x-coordinate and gamma_l

xlp = xlp_coord(2);

gamma_l = xlp - (1.0 - mu);

% normalize z-amplitude

az = amp_z / (dist * gamma_l);

c2 = cn_func(2, gamma_l);

c3 = cn_func(3, gamma_l);

c4 = cn_func(4, gamma_l);

lambda = sqrt((-c2 + 2.0 + sqrt(9.0 * c2^2 - 8.0 * c2)) / 2.0);

k = 2.0 * lambda/(lambda^2 + 1.0 - c2);

delta = lambda^2 - c2;

d1 = 16.0 * lambda^4 + 4.0 * lambda^2 * (c2 - 2.0) - 2.0 * c2^2 + c2 + 1.0;

d2 = 81.0 * lambda^4 + 9.0 * lambda^2 * (c2 - 2.0) - 2.0 * c2^2 + c2 + 1.0;

d3 = 2.0 * lambda * (lambda * (1.0 + k^2) - 2.0 * k);

a21 = 3.0 * c3 * (k^2 - 2.0) / 4.0 / (1.0 + 2.0 * c2);

a23 = -3.0 * lambda * c3 * (3.0 * k^3 * lambda - 6.0 * k * (k - lambda) ...
    + 4.0) / 4.0 / k / d1;

b21 = -3.0 * c3 * lambda * (3.0 * lambda * k - 4.0) / 2.0 / d1;

s1 = ((3.0 / 2.0) * c3 * (2.0 * a21 * (k^2 - 2.0) - a23 * (k^2 + 2.0)- ...
    2.0 * k * b21) - (3.0 / 8.0) * c4 * (3.0 * k^4 - 8.0 * k^2 + 8.0)) / d3;

a22 = 3.0 * c3 / 4.0 / (1.0 + 2.0 * c2);

a24 = -3.0 * c3 * lambda * (2.0 + 3.0 * lambda * k) / 4.0 / k / d1;

b22 = 3.0 * lambda * c3 / d1;

d21 = -c3 / 2.0 / lambda^2;

s2 = ((3.0 / 2.0) * c3 * (2.0 * a22 * (k^2 - 2.0) + ...
    a24 * (k^2 + 2.0) + 2.0 * k * b22 + 5.0 * d21) + ...
    (3.0 / 8.0) * c4 * (12.0 - k^2)) / d3;

a1 = -(3.0 / 2.0) * c3 * (2.0 * a21 + a23 + 5.0 * d21) - ...
    (3.0 / 8.0) * c4 * (12.0 - k^2);

a2 = (3.0 / 2.0) * c3 * (a24 - 2.0 * a22) + (9.0 / 8.0) * c4;

l1 = 2.0 * s1 * lambda^2 + a1;

l2 = 2.0 * s2 * lambda^2 + a2;

az_sqr = az^2;

% check if this az amplitude  is feasible

if (l1 ~= 0.0)
    
    term = (-l2 * az_sqr - delta) / l1;
    
else
    
    fprintf('\n\nerror in halo_sv - divide by zero\n\n');
    
    pause
    
end

ax = sqrt(term);

ax_sqr = ax^2;

% frequency correction

w = 1.0 + s1 * ax_sqr + s2 * az_sqr;

a31 = -9.0 * lambda * (c3 * (k * a23 - b21) + k * c4 * (1.0 + (1.0 / 4.0) ...
    * k^2)) / d2 + (9.0 * lambda^2 + 1.0 - c2) * (3.0 * c3 * (2.0 * a23 ...
    - k * b21) + c4 * (2.0 + 3.0 * k^2)) / 2.0 / d2;

a32 = -9.0 * lambda * (4.0 * c3 * (k * a24 - b22) + k * c4) / 4.0 / d2 - ...
    3.0 * (9.0 * lambda^2 + 1.0 - c2) * (c3 * (k * b22 + d21 - ...
    2.0 * a24) - c4) / 2.0 / d2;

b31 = (3.0 * lambda * (3.0 * c3 * (k * b21 - 2.0 * a23) - c4 * (2.0 + 3.0 * k^2)) + ...
    (9.0 * lambda^2 + 1.0 + 2.0 * c2) * (12.0 * c3 * (k * a23 - b21) + ...
    3.0 * k * c4 * (4.0 + k^2)) / 8.0) / d2;

b32 = (3.0 * lambda * (3.0 * c3 * (k * b22 + d21 - 2.0 * a24) - 3.0 * c4) + ...
    (9.0 * lambda^2 + 1.0 + 2.0 * c2) * (12.0 * c3 * (k * a24 - b22) + ...
    3.0 * c4 * k) / 8.0) / d2;

d31 = 3.0 * (4.0 * c3 * a24 + c4) / 64.0 / lambda^2;

d32 = 3.0 * (4.0 * c3 * (a23 - d21) + c4 * (4.0 + k^2)) / 64.0 / lambda^2;

% switch function

delta_n = 2 - hclass;

% normalized orbital period

period = 2.0 * pi / (lambda * w);

% compute normalized components of position vector

r_halo(1) = a21 * ax_sqr + a22 * az_sqr - ax * cos(t) + ...
    (a23 * ax_sqr - a24 * az_sqr) * cos(2.0 * t) + ...
    (a31 * ax^3 - a32 * ax * az_sqr) * cos(3.0 * t);

r_halo(2) = k * ax * sin(t) + (b21 * ax_sqr - b22 * az_sqr) * sin(2.0 * t) + ...
    (b31 * ax^3 - b32 * ax * az_sqr) * sin(3.0 * t);

r_halo(3) = delta_n * (az * cos(t) + d21 * ax * az * (cos(2.0 * t) - 3.0) + ...
    (d32 * az * ax_sqr - d31 * az^3) * cos(3.0 * t));

% compute normalized components of velocity vector

v_halo(1) = ax * sin(t) - (a23 * ax_sqr - a24 * az_sqr) * sin(2.0 * t) * 2.0 - ...
    (a31 * ax^3 - a32 * ax * az_sqr) * sin(3.0 * t) * 3.0;

v_halo(2) = k * ax * cos(t) + (b21 * ax_sqr - b22 * az_sqr) * cos(2.0 * t) * 2.0 + ...
    (b31 * ax^3 - b32 * ax * az_sqr) * cos(3.0 * t) * 3.0;

v_halo(3) = delta_n * (-az * sin(t) + d21 * ax * az * (-sin(2.0 * t) * 2.0)- ...
    (d32 * az * ax_sqr - d31 * az^3) * sin(3.0 * t) * 3.0);

% convert from richardson scale, libration point centered to
% standard normalized coordinates wrt barycenter

for i = 1:1:3
    
    r_crtbp(i) = r_halo(i) * gamma_l;
    
    v_crtbp(i) = v_halo(i) * gamma_l * (lambda * w);
    
end

r_crtbp(1) = r_crtbp(1) + xlp;



