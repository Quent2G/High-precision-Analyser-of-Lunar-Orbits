% HALO, High-precision Analyser of Lunar Orbits, is here built for the propagation of a lunar orbiter.
% The motion is integrated under 
% - the Keplerian attraction force
% - the perturbations due to the asymmetry of the lunar gravity field (up to 350*350 harmonics)
% - the perturbations due to the asymmetry of the Earth gravity field (up to 100*100 harmonics)
% - the third body attractions of the Sun and Jupiter
% - the solar radiation pressure
% - the General relativity
% - the Earth's albedo
%     
% The integration is performed in the "J2000" inertial frame, and the output is the state and time vectors.
%
% IMPORTANT: it is required to install mice toolbox from
% https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
%
% Authors:
% (c)Copyright 2017 - Ennio Condoleo
% Sapienza University of Rome, 08-25-2017
% ennio.condoleo@uniroma1.it
%
% 2024 - Quentin Granier
% Research project at University of New South Wales
% Supervisor's mail: yiyinfeixiong@gmail.com

    ts = cputime;

%% SPICE LIBRARIES
    addpath("mice/lib");
    addpath("mice/src/mice");
    addpath('prop');
    addpath('input');
    
    cspice_kclear; 
    metakernelcheck;
    cspice_furnsh('metakernel.tm');

%% FORCE MODEL
    orb = struct();
    orb = LoadModel("Ref",orb);
    
%% Initialization
    orb = LoadSequential(orb);
    Time = orb.seq.Time;

%%  Running the sequence
    SeqNames = fieldnames(orb.seq);
    options = odeset('RelTol',1e-7,'AbsTol',1e-14);    
    TStep = 30;

    for i = 2:length(SeqNames)
        Seq = orb.seq.(SeqNames{i});
        
        if Seq.type == "Propag"
            span = Seq.span;
            Tspan = Time:TStep:Time+span;
            [orb.seq.(SeqNames{i}).t,orb.seq.(SeqNames{i}).XJ2000]...
                = ode113(@prophpop,Tspan,orb.sat.X0iner,options,orb);

        elseif Seq.type == "DVPropag"
            X = orb.sat.X0iner;
            orb = LoadState(Seq.Orbi,orb);
            newX = orb.sat.X0iner;
            fprintf("Diff. pos. target = %f km\nDV2 = %f m/s\n",norm(newX(1:3)-X(1:3)),norm(newX(4:6)-X(4:6))*1e3)
            span = Seq.span;
            Tspan = Time:TStep:Time+span;
            [orb.seq.(SeqNames{i}).t,orb.seq.(SeqNames{i}).XJ2000]...
                = ode113(@prophpop,Tspan,[X(1:3) newX(4:6)],options,orb);

        elseif Seq.type == "fsolveProp"
            Period = Seq.T;
            Tspan = Time:TStep:Time+4*Period;
            X = fsolveConv(orb.sat.X0iner,orb);

            [orb.seq.(SeqNames{i}).t,orb.seq.(SeqNames{i}).XJ2000]...
                = ode113(@prophpop,Tspan,orb.sat.X0iner,options,orb);
            [~,orb.seq.(SeqNames{i}).XC] = ode113(@prophpop,Tspan,X,options,orb);

        elseif Seq.type == "Lambert"
            span = Seq.span;
            X0 = orb.sat.X0iner;
            r1 = X0(1:3);

            orb = LoadState(Seq.stop,orb);
            r2 = orb.sat.X0iner(1:3);

            [V1, V2, extremal_distances, exitflag] = lambert(r1, r2, span/86400, 0, orb.centralPlanet.GM);
            fprintf("DV1 = %f m/s\n",norm(V1-X0(4:6))*1e3)

            Tspan = Time:TStep:Time+Seq.span;
            [orb.seq.(SeqNames{i}).t,orb.seq.(SeqNames{i}).XJ2000]...
                = ode113(@prophpop,Tspan,[r1 V1],options,orb);
        
        elseif Seq.type == "OptimLambert"
            [t1,t2,X2,span] = LambertConv(orb); % Processing of the solution
            
            %------ propagation of the processed solution ------
            Tspan = Time:TStep:Time + t1;
            Time = Time+t1;
            [T1,XJ1] = ode113(@prophpop,Tspan,X0,options,orb);
            X1 = XJ1(end,:);

            [V1, V2, extremal_distances, exitflag] = lambert(X1(1:3), X2(1:3), span/86400, 0, orb.centralPlanet.GM);
            fprintf("DV1 = %f m/s\n",norm(V1-X1(4:6))*1e3) 
            Tspan = Time:TStep:Time + span;
            Time = Time+span;
            [T2,XJ2] = ode113(@prophpop,Tspan,[X1(1:3) V1],options,orb);
            X2t = XJ2(end,:);

            fprintf("Diff. pos. target = %f km\nDV2 = %f m/s\n",norm(X2t(1:3)-X2(1:3)),norm(V2-X2(4:6))*1e3)
            Tspan = Time:TStep:Time + t1;
            [T3,XJ3] = ode113(@prophpop,Tspan,[X2t(1:3) X2(4:6)],options,orb);
            
            span=t1; % to agree with next lines
            orb.seq.(SeqNames{i}).t = [T1;T2;T3];
            orb.seq.(SeqNames{i}).XJ2000 = [XJ1;XJ2;XJ3];
        end

        orb.sat.X0iner = orb.seq.(SeqNames{i}).XJ2000(end,:);
        Time = Time+span;
        orb.seq.Time = Time;
    end

    save('output/ORBdata','orb');
    disp("Run Time: "+num2str(cputime - ts)+"s")
    
%% Close Path
    rmpath('prop');
    rmpath('input');
    rmpath("mice/lib")
    rmpath("mice/src/mice")

