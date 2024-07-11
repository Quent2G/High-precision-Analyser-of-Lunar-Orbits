% set up options for ode45 for event finding

options = odeset('RelTol', 1.0e-12, 'AbsTol', 1.0e-12, 'Events', @halo_event);

tend = 100.0;

% solve for x-axis crossing (ry = 0)

[~, ~, tevent, yevent, ~] = ode45(@halo_eqms, [0.0 tend], [r_crtbp v_crtbp], options);