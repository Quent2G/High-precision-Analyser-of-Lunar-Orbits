function xlp = xlp_coord(ilp)

% location of libration point normalized x-component

% input

%  ilp = libration point of interest (1, 2 or 3)

% output

%  xlp = normalized x-coordinate of libration point

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convergence criterion for Brent's method

rtol = 1.0e-10;

% lower and upper bounds for coordinate search

xlwr = -2.0;
        
xupr = +2.0;
        
switch ilp
    
    case 1
        
        % L1 libration point (normalized)
        
        [xlp, ~] = brent('clp_func', xlwr, xupr, rtol);
        
    case 2
        
        % L2 libration point (normalized)
        
        [xlp, ~] = brent('clp_func', xlwr, xupr, rtol);
        
    case 3
        
        % L3 libration point (normalized)
        
        [xlp, ~] = brent('clp_func', xlwr, xupr, rtol);
        
end
