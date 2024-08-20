

function [t1,t2,X2,span] = LambertConv(orb)

    optionsfs = optimset('Display', 'iter', ...
        'MaxFunEval',100);
    
    % Tmin = [orb.seq.a.t1g(1) orb.seq.a.t2g(1) orb.seq.a.spang(1)];
    T0 = [orb.seq.a.t1g(2) orb.seq.a.t2g(2) orb.seq.a.spang(2)];
    % Tmax = [orb.seq.a.t1g(3) orb.seq.a.t2g(3) orb.seq.a.spang(3)];
    T = fminsearch(@(T)DVTOT(T,orb),T0,optionsfs);
    [t1,t2,span] = deal(T(1),T(2),T(3));

    target = orb.seq.a.target;
    orb2 = LoadState(target,orb);
    X2 = Pos(t2,orb2.sat.X0iner,orb);
end

function DV = DVTOT(T,orb)
    X1 = Pos(T(1),orb.sat.X0iner,orb);
    
    target = orb.seq.a.target;
    orb2 = LoadState(target,orb);
    X2 = Pos(T(2),orb2.sat.X0iner,orb);

    [V1, V2, ~, ~] = lambert(X1(1:3), X2(1:3), T(3)/86400, 0, orb.centralPlanet.GM);
    
    %Assessment good Lambert
    X2t = Pos(T(3),[X1(1:3) V1],orb);
    pun = norm(X2t(1:3)-X2(1:3))*1/2e3; %Punishment
    % pun=0;
    
    DV1 = norm(X1(4:6)-V1);
    DV2 = norm(X2(4:6)-V2);
    DV = DV1 + DV2 + pun;
end

function X = Pos(t,X0,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    [~,Y] = ode45(@prophpop,[0 t],X0,options,orb);
    X = Y(end,:);
end