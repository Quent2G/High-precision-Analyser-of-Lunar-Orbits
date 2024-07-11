function y = clp_func(x)

% normalized x-coordinate function

% input

%  x = current argument

% output

%  y = objective function evaluated at x

% global

%  ilp = libration point number (1, 2 or 3)
%  mu  = normalized gravity constant

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ilp mu

switch ilp
    
    case 1
        
        % L1 libration point
        
        y = x - (1.0 - mu) / (mu + x)^2 + mu / (x - 1.0 + mu)^2;
        
    case 2
        
        % L2 libration point
        
        y = x - (1.0 - mu) / (mu + x)^2 - mu / (x - 1.0 + mu)^2;
        
    case 3
        
        % L3 libration point
        
        y = x + (1.0 - mu) / (mu + x)^2 + mu / (x - 1.0 + mu)^2;
        
end
