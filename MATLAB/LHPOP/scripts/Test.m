% Define initial guess
r0 = [7000; 0]; % Initial position (km)
v0 = [0; 7.5];  % Initial velocity (km/s)
initial_state = [r0; v0];

% Define target final state
r_target = [42164; 0]; % Target position (km)
v_target = [0; 3.07];  % Target velocity (km/s)
target_state = [r_target; v_target];

% Define tolerance and maximum iterations
tolerance = 1;
max_iterations = 100;

% Time of flight
tof = 3600; % 1 hour (s)

% Function to propagate state using simple two-body dynamics
propagate_state = @(state, dt) [state(1:2) + state(3:4) * dt; state(3:4)];

% Initial correction
correction = zeros(size(initial_state));

for iteration = 1:max_iterations
    % Apply correction to initial state
    corrected_state = initial_state + correction;
    disp(corrected_state);
    
    % Propagate corrected state to final time
    final_state = propagate_state(corrected_state, tof);
    
    % Compute error in final state
    error = target_state - final_state;
    
    % Check convergence
    if norm(error) < tolerance
        fprintf('Converged in %d iterations.\n', iteration);
        break;
    end
    
    % Compute Jacobian (partial derivatives of final state with respect to initial state)
    dt = 1e-6;
    jacobian = zeros(length(final_state), length(initial_state));
    for i = 1:length(initial_state)
        perturbation = zeros(size(initial_state));
        perturbation(i) = dt;
        perturbed_final_state = propagate_state(corrected_state + perturbation, tof);
        jacobian(:, i) = (perturbed_final_state - final_state) / dt;
    end
    
    % Compute correction using minimum-norm solution
    correction = correction + jacobian' * (jacobian * jacobian')^(-1) * error;
end

% Check if the solution converged
if norm(error) >= tolerance
    fprintf('Did not converge within the maximum number of iterations.\n');
end

% Display final corrected state and error
final_state = propagate_state(initial_state + correction, tof);
error = target_state - final_state;
fprintf('Final corrected state:\n');
disp(final_state);
fprintf('Final error:\n');
disp(error);
