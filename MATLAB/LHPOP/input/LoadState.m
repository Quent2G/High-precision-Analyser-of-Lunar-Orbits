function [orb] = LoadState(Orbit,orb)

    if Orbit == "NRHO"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
        
        X = -4103.547489327654; % Semi Major Axis
        Y = 20820.76103083731;    % Eccentricity
        Z = -68714.1843688717;    % Inclination
        VX = -0.0641922245925115;    % Longitude of Ascending Node
        VY = -0.004286985491924689;    % Argument of Periapsis
        VZ  = 0.0022994082206247227;    % Mean Anomaly

    elseif Orbit == "DRO14"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
        
        X = -1642.1842792591378; % Semi Major Axis
        Y = 62149.566213382524;    % Eccentricity
        Z = 32159.59842593645;    % Inclination
        VX = 0.31558332872434935;    % Longitude of Ascending Node
        VY = 0.015087340226611998;    % Argument of Periapsis
        VZ  = -0.014392035576109786;    % Mean Anomaly
    
    elseif Orbit == "DRO14C"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
        
        X = -25423.19272080875; % Semi Major Axis
        Y = 63798.47179482236;    % Eccentricity
        Z = 34566.169617323954;    % Inclination
        VX = 0.3078236588788274;    % Longitude of Ascending Node
        VY = -0.04046560323071425;    % Argument of Periapsis
        VZ  = -0.04225114303350116;    % Mean Anomaly

    elseif Orbit == "ELFO"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 12.931; % Semi Major Axis
        Y = 235.558;    % Eccentricity
        Z = 2388.38;    % Inclination
        VX = -1.80519;    % Longitude of Ascending Node
        VY = 0.0991;    % Argument of Periapsis
        VZ  = 3.83e-07;    % Mean Anomaly

    elseif Orbit == "CLEM15"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -2098.4745703999997; % Semi Major Axis
        Y = 1298.5243538999998;    % Eccentricity
        Z = -3494.2586871999997;    % Inclination
        VX = 0.8948704814799997;    % Longitude of Ascending Node
        VY = 0.17693000307000004;    % Argument of Periapsis
        VZ  = -0.15601097269;    % Mean Anomaly

    elseif Orbit == "Capstone25Nov"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -16983.14075642353; % Semi Major Axis
        Y = 21213.5584242304;    % Eccentricity
        Z = -58035.6304537942;    % Inclination
        VX = -0.04457645856905286;    % Longitude of Ascending Node
        VY = -0.05137482728591616;    % Argument of Periapsis
        VZ  = 0.1416421540429468;    % Mean Anomaly
    
    elseif Orbit == "Orion26Nov"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 91935.04138932908; % Semi Major Axis
        Y = -11116.52719922882;    % Eccentricity
        Z = -12261.72337590001;    % Inclination
        VX = -0.08994654608636954;    % Longitude of Ascending Node
        VY = -0.08923662419613348;    % Argument of Periapsis
        VZ  = -0.03967967460594898;    % Mean Anomaly

    elseif Orbit == "Orion29Nov"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 41570.98496247831; % Semi Major Axis
        Y = -50015.47166299814;    % Eccentricity
        Z = -28760.46002667214;    % Inclination
        VX = -0.2465381485425481;    % Longitude of Ascending Node
        VY = -0.182296424457602;    % Argument of Periapsis
        VZ  = -0.07642694222326024;    % Mean Anomaly

    elseif Orbit == "Orion29Nov16"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 27710.37806303968; % Semi Major Axis
        Y = -60030.4657139458;    % Eccentricity
        Z = -32932.94443647428;    % Inclination
        VX = -0.2314339589046175;    % Longitude of Ascending Node
        VY = -0.1626823717943428;    % Argument of Periapsis
        VZ  = -0.06728693114754408;    % Mean Anomaly

    elseif Orbit == "LRO6Feb"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -514.825148; % Semi Major Axis
        Y = -183.032603;    % Eccentricity
        Z = -1734.88976;    % Inclination
        VX = -0.595944951;    % Longitude of Ascending Node
        VY = -1.49949921;    % Argument of Periapsis
        VZ  = 0.331244402;    % Mean Anomaly
    
    elseif Orbit == "LRO1Feb"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 820.1514640535461; % Semi Major Axis
        Y = 1091.2367697664858;    % Eccentricity
        Z = 1262.7894519284164;    % Inclination
        VX = 0.22458861495091112;    % Longitude of Ascending Node
        VY = 1.1177228435015545;    % Argument of Periapsis
        VZ  = -1.1391038535297986;    % Mean Anomaly
    
    elseif Orbit == "Perso"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -175.001861; % Semi Major Axis
        Y = -1175.17121;    % Eccentricity
        Z = 1410.69401;    % Inclination
        VX = 0.791300523;    % Longitude of Ascending Node
        VY = 1.03420255;    % Argument of Periapsis
        VZ  = 0.971046194;    % Mean Anomaly
    end

%% Initialization
% position and velocity
    % =========================================================================% 
    orb.frame.initstate = RefS;
    if strcmp(CoordT,'Cartesian')
        orb.sat.X0 = [X,Y,Z,VX,VY,VZ];
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_sxform(orb.frame.initstate,orb.frame.to,orb.epoch.et(1));
            orb.sat.X0iner = Rgeog_iner*orb.sat.X0;
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    else
        orb.sat.keplstate(1,1:2) = {'SMA',SMA};        % semi-major axis             [km]
        orb.sat.keplstate(2,1:2) = {'ECC',ECC};        % eccentricity
        orb.sat.keplstate(3,1:2) = {'INC',INC*(pi/180)};  % inclination                 [rad]
        orb.sat.keplstate(4,1:2) = {'LAN',LAN*(pi/180)};   % longitude of ascending node [rad]
        orb.sat.keplstate(5,1:2) = {'AOP',AOP*(pi/180)};   % argument of pericenter      [rad]
        orb.sat.keplstate(6,1:2)  = {'MA',MA*(pi/180)};
        orb.sat.X0 = (cspice_conics([ orb.sat.keplstate{1,2}*(1-orb.sat.keplstate{2,2}),orb.sat.keplstate{2,2},orb.sat.keplstate{3,2},...
                                    orb.sat.keplstate{4,2},orb.sat.keplstate{5,2},orb.sat.keplstate{6,2},...
                                    orb.epoch.et(1),orb.centralPlanet.GM]',orb.epoch.et(1)))';
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_pxform(orb.frame.initstate,orb.frame.to,orb.epoch.et(1));
            orb.sat.X0iner(1:3,1) = Rgeog_iner*orb.sat.X0(1,1:3)';
            orb.sat.X0iner(4:6,1) = Rgeog_iner*orb.sat.X0(1,4:6)';
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    end
    % =========================================================================% 

end

