

function X = fsolveConv(Xi,orb)
    FDSS=[1e4,0.1];
    optionsfs = optimoptions('fsolve','Display', 'iter', ...
        'MaxFunctionEvaluations',200);
        % 'FiniteDifferenceStepSize',[FDSS(1),FDSS(1),FDSS(1),FDSS(2),FDSS(2),FDSS(2)], ...
    
    for i = 1:1
        X = fsolve(@(X)Propag(X,orb),Xi,optionsfs);
        % options = odeset('RelTol',1e-7,'AbsTol',1e-12);
        % [~,Ytraj] = ode45(@prophpop,orb.epoch.et,X,options,orb);
        % Xi=Ytraj(end,:);
    end
end

function [Mri,Mir] = MatRI(t)
    SE = cspice_spkezr('EARTH',t,'J2000','NONE','MOON');

    uR = -SE(1:3)/norm(SE(1:3));
    uTh = -SE(4:6)/norm(SE(4:6));
    uZ = cross(uR,uTh);

    Mri = [uR,uTh,uZ];
    Mir = inv(Mri);
end

function dY = Propag(X,orb)
    options = odeset('RelTol',1e-7,'AbsTol',1e-12);
    span = orb.epoch.et+[0,orb.epoch.et(2)-orb.epoch.et(1)];
    [~,Ytraj] = ode45(@prophpop,span,X,options,orb);
    Y = Ytraj(end,:);
    [~,Miri] = MatRI(span(1));
    [~,Mirf] = MatRI(span(2));
    Xr = Miri*X(1:3)';
    Yr = Mirf*Y(1:3)';
    Xpr = Miri*X(4:6)';
    Ypr = Mirf*Y(4:6)';
    dY = [Yr' - Xr',Ypr' - Xpr'];
end