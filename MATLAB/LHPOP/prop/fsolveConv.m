

function X = fsolveConv(Xi,orb)
    % FDSS=[1e4,0.1];
    optionsfs = optimoptions('fsolve','Display', 'iter', ...
        'MaxFunctionEvaluations',200);
        % 'FiniteDifferenceStepSize',[FDSS(1),FDSS(1),FDSS(1),FDSS(2),FDSS(2),FDSS(2)], ...
    
    X = fsolve(@(X)Propag(X,orb),Xi,optionsfs);
    % V = fsolve(@(V)Propag2(V,orb),[Xi orb.seq.Time+orb.seq.a.span Xi orb.seq.Time+2*orb.seq.a.span],optionsfs);
    % X = V(1:6);
end

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
    Omg = (uR-uRm)/dt/uTh;

    Xr(1:3) = Mir*X(1:3)';
    Xr(4:6) = Mir*(X(4:6)' - cross(Omg*uZ,X(1:3)'));
end

function dY = Propag(X,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    n=3;                    %number of periods taken into account

    % span = orb.epoch.et+[0,orb.epoch.et(2)-orb.epoch.et(1)];
    span = zeros(1,n+1);
    for i = 0:n
        span(i+1) = orb.seq.Time+i*orb.seq.a.span;
    end

    [~,Ytraj] = ode45(@prophpop,span,X,options,orb);

    [~,~,Xr] = MatRI(span(1),X);

    dY = zeros(6*n,1);
    for i = 0:n-1
        [~,~,Ytr] = MatRI(span(i+2),Ytraj(i+2,:));
        dY(6*i + 1:6*i + 3) = Ytr(1:3)'-Xr(1:3)';
        dY(6*i + 4:6*i + 6) = Ytr(4:6)'-Xr(4:6)';
    end
end


function dY = Propag2(V,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    dY = zeros(12);

    t0 = orb.seq.Time;
    S0 = V(1:6);
    t1 = V(7);
    S1 = V(8:13);
    t2 = V(14);

    span = [t0 t1];
    [~,Ytraj] = ode45(@prophpop,span,S0,options,orb);
    P0 = Ytraj(2,:);
    [~,~,S1r] = MatRI(t1,S1);
    [~,~,P0r] = MatRI(t1,P0);
    dY(1:6) = P0r-S1r;

    span = [t1 t2];
    [~,Ytraj] = ode45(@prophpop,span,S1,options,orb);
    P1 = Ytraj(2,:);
    [~,~,P1r] = MatRI(t1,P1);
    [~,~,S0r] = MatRI(t0,S0);
    dY(1:6) = P1r-S0r;
end