function [Xi] = DiffCorrector(m,prop,orb,options)
%DIFFCORRECTOR Summary of this function goes here
%   Detailed explanation goes here
    
    Xi = orb.sat.X0iner;
    for k = 1:m
        % Jacobian computation
        n=6;
        J = zeros(n, n);
        I = eye(6);
        h = [1e-3,1e-6];             %order m and mm/s
        [~,Y] = ode45(@prophpop,orb.epoch.span,Xi,options,orb);
        Yref = Y(end,:);
    
        for i = 1:n
            step = h(1+(i>3));
            X = Xi + I(i,:)*step;
            [~,Y] = ode45(prop,orb.epoch.span,X,options,orb);
            J(:, i) = (Y(end,:) - Yref) / step;
        end

        % dX = J\(orb.XJ2000(1,:)-Yref)';
        dX = J'*(J*J')^(-1)*(orb.XJ2000(1,:)-Yref)';
        Xi = Xi + dX';
    end
end

