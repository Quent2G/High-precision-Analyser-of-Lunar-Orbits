def ConvMatlab(PosInit):
    print(f"""X = {PosInit[0]}; % Semi Major Axis
        Y = {PosInit[1]};    % Eccentricity
        Z = {PosInit[2]};    % Inclination
        VX = {PosInit[3]};    % Longitude of Ascending Node
        VY = {PosInit[4]};    % Argument of Periapsis
        VZ  = {PosInit[5]};    % Mean Anomaly""")