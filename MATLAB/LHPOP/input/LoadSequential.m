function [orb] = LoadSequential(orb)
%LOADSEQUENTIAL Summary of this function goes here
%   Detailed explanation goes here

    orb = LoadState("Orion29Nov16",orb);
    Time = '2022 Nov 29 16:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    % Time = '2020 Feb 06 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    % Time = '1994 Apr 15 15:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    orb.seq.Time = cspice_str2et(Time);

    % orb.seq.a.type = "fsolveProp";
    orb.seq.a.type = "Propag";
    orb.seq.a.span = 1*86400;
    % orb.seq.a.span = 5.7444e+05;
    % orb.seq.a.span = 5.7444e+05*0.45;
    % orb.seq.a.span = 3600*15;
    % orb.seq.a.span = 5.7444e+05*1.1;

    % orb.seq.a.type = "OptimLambert";
    % orb.seq.a.target = "ELFO";
    % orb.seq.a.t1g = [5.7444e+05, 5.7444e+05*1.1, 5.7444e+05*2];
    % orb.seq.a.t2g = [3600*6 3600*11 3600*18];
    % orb.seq.a.spang = [3600 3600*6 2*86400];

    % orb.seq.b.type = "Lambert";
    % orb.seq.b.stop = "ELFO";
    % orb.seq.b.span = 3600*11;
    % 
    % orb.seq.c.type = "DVPropag";
    % orb.seq.c.Orbi = "ELFO";
    % orb.seq.c.span = 3600*6;
end

