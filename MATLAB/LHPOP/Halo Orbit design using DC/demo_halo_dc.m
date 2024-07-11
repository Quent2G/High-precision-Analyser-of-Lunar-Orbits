% demo_halo_dc.m      February 21, 2022

% demonstrates halo orbit design in the Earth-Moon
% system using differential correction (dc)

% class I halo orbit about the L2 libration point

% script features

% (1) initial guess using Richardson's algorithm
% (2) solve system of nonlinear equations
% (3) graphics display of dc-refined halo orbit

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu hclass dist yi tevent yevent

y = zeros(6, 1);

clc; home;

fprintf('\ndemo_halo_dc - halo orbit design using differential correction\n\n');

while(1)

    fprintf('\nplease input the z-amplitude in kilometers\n');

    ampl_z = input('? ');

    if (ampl_z > 0.0)

        break;

    end

end

% create z-amplitude string for display

az_str = num2str(ampl_z, '(z amplitude = %8.2f kilometers)\n');

% create class I trajectories

hclass = 1;

% -----------------
% earth-moon system
% -----------------

% gravitational constant of earth (kilometers^3/seconds^2)

mu1 = 398600.4415;

% gravitational constant of the moon (kilometers^3/seconds^2)

mu2 = 4902.8;

% normalized gravitational constant

mu = mu2 / (mu1 + mu2);

% distance between primary and secondary bodies (kilometers)

dist = 384400.0;

% length unit (kilometers)

lu = dist;

% mean motion (radians/second)

mmotion = sqrt((mu1 + mu2) / lu^3);

% time unit (seconds)

tu = 1.0 / mmotion;

% compute normalized coordinates of L2 libration point

[xl2, yl2] = lp_coordinates;

% -------------------------------------------------
% compute initial guess for normalized state vector
% and orbital period using richardson's algorithm
% -------------------------------------------------

[r_halo, v_halo, period] = halo_sv(ampl_z, 0.0);

for i = 1:1:3
    
    yi(i) = r_halo(i);
    
    yi(i + 3) = v_halo(i);
    
end

% -----------------------------------
% solve system of nonlinear equations
% -----------------------------------

% initial guess for normalized x-component of position

xg(1) = v_halo(1);

% initial guess for normalized y-component of velocity

xg(2) = v_halo(3);

% define algorithm options

options = optimoptions('fsolve', 'algorithm', 'trust-region-dogleg');

% fsolve

[xsol, fofx, exitflag] = fsolve(@halo_objective, xg, options);

% display initial guess from richardson's algorithm

fprintf('\ninitial normalized state vector - richardson algorithm');
fprintf('\n------------------------------------------------------\n');

svprint(yi(1:3), yi(4:6));

% ------------------------------
% set ode45 options for plotting
% ------------------------------

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);

% initialize

dt = 0.01;

t2 = -dt;

npts = fix(2.0 * tevent / dt) + 1;

xplot = zeros(npts, 1);

yplot = zeros(npts, 1);

zplot = zeros(npts, 1);

% dc-corrected normalized state vector

for i = 1:1:6
    
    y(i) = yi(i);
    
end

y(4) = xsol(1);

y(6) = xsol(2);

y_dc = y;

% display dc-refined initial conditions

fprintf('\ndc-refined initial normalized state vector');
fprintf('\n------------------------------------------\n');

svprint(y(1:3), y(4:6));

% display dc-refined time and state vector at y-axis crossing

fprintf('\ndc-refined normalized state vector at y-axis crossing');
fprintf('\n-----------------------------------------------------\n');

fprintf('\ntevent (days) =  %14.10f\n', tevent * tu / 86400.0);

svprint(yevent(1:3), yevent(4:6));

% initial dimensional position data for plotting

xplot(1) = dist * y(1);

yplot(1) = dist * y(2);

zplot(1) = dist * y(3);

% compute integrated dimensional position data for plotting

for i = 2:1:npts
    
    t1 = t2;
    
    t2 = t1 + dt;
    
    [twrk, ysol] = ode45(@halo_eqms, [t1, t2], y, options);
    
    xplot(i) = dist * ysol(length(twrk), 1);
    
    yplot(i) = dist * ysol(length(twrk), 2);
    
    zplot(i) = dist * ysol(length(twrk), 3);
    
    y = ysol(length(twrk), 1:6);
    
end

% display time and state vector after intergating the 
% dc-refined initial conditions for one orbital period

fprintf('\ndc-refined normalized state vector integrated for one period');
fprintf('\n------------------------------------------------------------\n');

fprintf('\ntevent (days) =  %14.10f\n', 2.0 * tevent * tu / 86400.0);

svprint(y(1:3), y(4:6));

% --------------------------------------
% create Halo orbit display in x-y plane
% --------------------------------------

figure(1);

hold on;

grid on;

axis equal;

plot(xplot, yplot, '-k', 'LineWidth', 1.5);

plot(xplot(1), yplot(1), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

xlabel('x coordinate (kilometers)', 'FontSize', 14);

ylabel('y coordinate (kilometers)', 'FontSize', 14);

% plot location of moon (blue)

plot(dist * (1.0 - mu), 0.0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% plot location of L2 libration point (red)

plot(dist * xl2, dist * yl2, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

title({'L_2 Halo Orbit Geometry', az_str},'FontSize', 16);

% create tiff graphics disk file

print ('halo_dc_xy.tif', '-dtiff');

% --------------------------------------
% create Halo orbit display in x-z plane
% --------------------------------------

figure(2);

hold on;

grid on;

axis equal;

plot(xplot, zplot, '-k', 'LineWidth', 1.5);

plot(xplot(1), zplot(1), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

xlabel('x coordinate (kilometers)', 'FontSize', 14);

ylabel('z coordinate (kilometers)', 'FontSize', 14);

% plot location of moon (blue)

plot(dist * (1.0 - mu), 0.0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% plot location of L2 libration point (red)

plot(dist * xl2, dist * yl2, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

title({'L_2 Halo Orbit Geometry', az_str},'FontSize', 16);

% create tiff graphics disk file

print ('halo_dc_xz.tif', '-dtiff');

% --------------------------------------
% create Halo orbit display in y-z plane
% --------------------------------------

figure(3);

hold on;

grid on;

axis equal;

plot(yplot, zplot, '-k', 'LineWidth', 1.5);

plot(yplot(1), zplot(1), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

xlabel('y coordinate (kilometers)', 'FontSize', 14);

ylabel('z coordinate (kilometers)', 'FontSize', 14);

% plot location of moon (blue)

plot(0.0, 0.0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% plot location of L2 libration point (red)

plot(dist * yl2, 0.0, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

title({'L_2 Halo Orbit Geometry', az_str},'FontSize', 16);

% create tiff graphics disk file

print ('halo_dc_yz.tif', '-dtiff');

% -------------------------------------------
% create halo orbit three-dimensional display
% -------------------------------------------

figure(4);

hold on;

grid on;

axis equal;

plot3(xplot, yplot, zplot, '-k', 'LineWidth', 1.5);

plot3(xplot(1), yplot(1), zplot(1), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

xlabel('x coordinate (kilometers)', 'FontSize', 14);

ylabel('y coordinate (kilometers)', 'FontSize', 14);

zlabel('z coordinate (kilometers)', 'FontSize', 14);

% plot location of the moon (blue)

plot3(dist * (1.0 - mu), 0.0, 0.0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% plot location of L2 libration point (red)

plot3(dist * xl2, 0.0, 0.0, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

title({'L_2 Halo Orbit Geometry', az_str},'FontSize', 16);

rotate3d on;

% set view

view(-15.0, 50.0);

% create tiff graphics disk file

print('halo_dc_3d.tif', '-dtiff');

fprintf('\n\n');


