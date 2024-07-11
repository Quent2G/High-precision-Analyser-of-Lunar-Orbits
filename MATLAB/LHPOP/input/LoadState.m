function [orb] = LoadState(Orbit,orb)

    if Orbit == "NRHO"
        X = -4103.547489327654; % Semi Major Axis
        Y = 20820.76103083731;    % Eccentricity
        Z = -68714.1843688717;    % Inclination
        VX = -0.0641922245925115;    % Longitude of Ascending Node
        VY = -0.004286985491924689;    % Argument of Periapsis
        VZ  = 0.0022994082206247227;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2022 Dec 2 15:34:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

    elseif Orbit == "NRHOl"
        X = -4103.547489327654; % Semi Major Axis
        Y = 20820.76103083731;    % Eccentricity
        Z = -68714.1843688717;    % Inclination
        VX = -0.0641922245925115;    % Longitude of Ascending Node
        VY = -0.004286985491924689;    % Argument of Periapsis
        VZ  = 0.0022994082206247227;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2022 Dec 26 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

    elseif Orbit == "DRO14"
        X = -1642.1842792591378; % Semi Major Axis
        Y = 62149.566213382524;    % Eccentricity
        Z = 32159.59842593645;    % Inclination
        VX = 0.3549830252637374;    % Longitude of Ascending Node
        VY = 0.017074626095640213;    % Argument of Periapsis
        VZ  = -0.016135387025700578;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2022 Dec 24 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

    elseif Orbit == "DRO14l"
        X = -1642.1842792591378; % Semi Major Axis
        Y = 62149.566213382524;    % Eccentricity
        Z = 32159.59842593645;    % Inclination
        VX = 0.3549830252637374;    % Longitude of Ascending Node
        VY = 0.017074626095640213;    % Argument of Periapsis
        VZ  = -0.016135387025700578;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2023 Mar 26 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

    elseif Orbit == "NRHOConv"
        X = -4420.834072073623; % Semi Major Axis
        Y = 20565.126095586234;    % Eccentricity
        Z = -68680.07125961124;    % Inclination
        VX = -0.05260322792429791;    % Longitude of Ascending Node
        VY = 0.005711627096215483;    % Argument of Periapsis
        VZ  = 0.012345164825294238;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2022 Dec 21 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
    
    elseif Orbit == "Perso"
        X = -4420.834072073623; % Semi Major Axis
        Y = 20565.126095586234;    % Eccentricity
        Z = -68680.07125961124;    % Inclination
        VX = -0.05260322792429791;    % Longitude of Ascending Node
        VY = 0.005711627096215483;    % Argument of Periapsis
        VZ  = 0.012345164825294238;    % Mean Anomaly

        TStart = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
        TStop =  '2022 Dec 21 00:00:00.000';  % stop time [yyyy month dd hh:mm:ss.---]
        TStep = 30; % Time step [s]
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
    end

%% Analisys start time and stop time
    % ==================================================================================================================================================
    orb.epoch.start = TStart;
    orb.epoch.stop  = TStop;
    
    orb.epoch.str   = char(orb.epoch.start,orb.epoch.stop);   
    orb.epoch.et    = cspice_str2et(orb.epoch.str); % seconds past the Julian Date
    
    orb.epoch.span = orb.epoch.et(1):TStep:orb.epoch.et(end);%+params.orb.periodTrig;
    % ==================================================================================================================================================

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

