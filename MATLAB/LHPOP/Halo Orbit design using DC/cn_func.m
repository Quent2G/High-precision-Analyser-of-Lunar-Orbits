function cn = cn_func(n, gl)

% halo_sv c_n support function

% input

%  n = c_n identifier

% output

%  cn = c_n function evaluated for argument n

% global

%  ilp = libration point number
%  mu  = normalized gravity constant

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ilp mu

switch ilp
    
    case 1
        
        % L1 libration point
        
        cn = (mu + (-1.0)^n * (1.0 - mu) * (gl / (1.0 - gl))^(n + 1)) / gl^3;
        
    case 2
        
        % L2 libration point
        
        cn = ((-1)^n * (mu + (1.0 - mu) * (gl / (1.0 + gl))^(n + 1))) / gl^3;
        
    case 3
        
        % L3 libration point
        
        cn = (1.0 - mu + mu * (gl / (1.0 + gl))^(n + 1)) / gl^3;
        
end
