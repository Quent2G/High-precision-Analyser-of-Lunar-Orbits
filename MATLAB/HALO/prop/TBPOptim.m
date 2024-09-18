% TBPOptim allows to optimize the initial state of a 3-body problem orbit in
% order to get an ephemeris 3BP orbit from a good initial guess, like a
% CR3BP solution fitted to ephemeris.

function X = TBPOptim(Xi,orb)
    optionsfs = optimoptions('fsolve','Display', 'iter', ...
        'MaxFunctionEvaluations',200);
    
    % fsolve tries to nullify Propag
    X = fsolve(@(X)Propag(X,orb),Xi,optionsfs);
end

% Propag returns a list of 6*n equations to nullify, n being the number of periods taken into account
% For each period T, the 6 equations are the difference of each component of
% the initial state and the propagated state at period T.
function dY = Propag(X,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-14);
    n=3;

    span = zeros(1,n+1);
    for i = 0:n
        span(i+1) = orb.seq.Time+i*orb.seq.a.T;
    end

    [~,Ytraj] = ode113(@prophpop,span,X,options,orb);

    [~,~,Xr] = MatRI(span(1),X);

    dY = zeros(6*n,1);
    for i = 0:n-1
        [~,~,Ytr] = MatRI(span(i+2),Ytraj(i+2,:));
        dY(6*i + 1:6*i + 3) = Ytr(1:3)'-Xr(1:3)';
        dY(6*i + 4:6*i + 6) = Ytr(4:6)'-Xr(4:6)';
    end
end

% MatRI computes the transition matrix from rotational to inertial at a given time, its
% inverse and the vector Xr which is X in the rotational frame.
function [Mri,Mir,Xr] = MatRI(t,X)
    SE = cspice_spkezr('EARTH',t,'J2000','NONE','MOON');
    dt = 1;
    SEm = cspice_spkezr('EARTH',t-dt,'J2000','NONE','MOON');

    uR = -SE(1:3)/norm(SE(1:3));
    uRm = -SEm(1:3)/norm(SEm(1:3));
    uTh = -SE(4:6)/norm(SE(4:6));
    uZ = cross(uR,uTh);

    Mri = [uR,uTh,uZ];
    Mir = inv(Mri);
    Omg = norm(uR-uRm)/dt/norm(uTh);

    Xr(1:3) = Mir*X(1:3)';
    Xr(4:6) = Mir*(X(4:6)' - cross(Omg*uZ,X(1:3)'));
end