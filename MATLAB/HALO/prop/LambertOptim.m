% LambertOptim finds the best lambert transfert between two orbits.
% The mission will then be a free propagation on the first
% orbit until best start, then lambert and change of
% orbit, and then free propagation of same length as the
% first one.

function [t1,t2,X2,span] = LambertOptim(orb)

    optionsfs = optimset('Display', 'iter', ...
        'MaxFunEval',100);
    
    % fminsearch tries to minimise DVTOT
    T0 = [orb.seq.a.t1 orb.seq.a.t2 orb.seq.a.span];
    T = fminsearch(@(T)DVTOT(T,orb),T0,optionsfs);
    [t1,t2,span] = deal(T(1),T(2),T(3));

    target = orb.seq.a.target;
    orb2 = LoadState(target,orb);
    X2 = Pos(t2,orb2.sat.X0iner,orb);
end

% DVTOT gives the total DV of the lambert transfert modified by a
% punishment term of the error to target is bad.
function DV = DVTOT(T,orb)
    X1 = Pos(T(1),orb.sat.X0iner,orb);
    
    target = orb.seq.a.target;
    orb2 = LoadState(target,orb);
    X2 = Pos(T(2),orb2.sat.X0iner,orb);

    [V1, V2, ~, ~] = lambert(X1(1:3), X2(1:3), T(3)/86400, 0, orb.centralPlanet.GM);
    
    %Assessment good Lambert
    X2t = Pos(T(3),[X1(1:3) V1],orb);
    pun = norm(X2t(1:3)-X2(1:3))*1/2e3; %Punishment
    
    DV1 = norm(X1(4:6)-V1);
    DV2 = norm(X2(4:6)-V2);
    DV = DV1 + DV2 + pun;
end

% Pos propagates an initial state X0 for a time t
function X = Pos(t,X0,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    [~,Y] = ode113(@prophpop,[0 t],X0,options,orb);
    X = Y(end,:);
end