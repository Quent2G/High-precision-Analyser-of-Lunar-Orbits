function [X] = prophpopGraphAccT(t,X0,model)
    
    global TextAcc;
    fprintf(TextAcc, "t   "+num2str(t)+"\n");

    % state in J2000 reference system centred at the central planet
    xJ2000 = X0(1:3);
    vJ2000 = X0(4:6);
    
    % acceleration due to the attraction of other planets
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_pointMasses = accelpntmassesGraph(xJ2000,model.pointMasses.stringName,model.pointMasses.GM,...    
        t,model.frame.integr,model.centralPlanet.stringName,model);
    % ---------------------------------------------------------------------------------------------------------- %

    % acceleration due to the central planet gravity field
    % ---------------------------------------------------------------------------------------------------------- %
    Rgeog_iner = cspice_pxfrm2(model.frame.from,model.frame.to,t,t);
    acc_centralPlanet = accelharmonic(xJ2000,Rgeog_iner',model.prop.harmonics.degree,model.prop.harmonics.order,...
        model.prop.harmonics.Cnm,model.prop.harmonics.Snm,model.centralPlanet.GM,model.centralPlanet.RE);
    fprintf(TextAcc, "LunGrav   "+num2str(norm(acc_centralPlanet))+"\n");
    % ---------------------------------------------------------------------------------------------------------- %   
    
    % acceleration due to the Earth gravity field
    % ---------------------------------------------------------------------------------------------------------- %
    REarth_iner = cspice_pxfrm2(model.frame.fromE,model.frame.to,t,t);
    R_Moon_Earth = cspice_spkezr('EARTH',t,model.frame.to,'NONE','MOON');
    R_Moon_Earth = R_Moon_Earth(1:3);
    acc_earth = accelharmonic(xJ2000-R_Moon_Earth,REarth_iner',model.prop.harmonics.degreeE,model.prop.harmonics.orderE,...
                model.prop.harmonics.ECnm,model.prop.harmonics.ESnm,model.Earth.GM,model.Earth.RE) + ...
                -accelharmonic(-R_Moon_Earth,REarth_iner',model.prop.harmonics.degreeE,model.prop.harmonics.orderE,...
                model.prop.harmonics.ECnm,model.prop.harmonics.ESnm,model.Earth.GM,model.Earth.RE);
    fprintf(TextAcc, "EarthGrav   "+num2str(norm(acc_earth))+"\n");
    % ---------------------------------------------------------------------------------------------------------- %   
    
    % acceleration due to the general relativity
    % ---------------------------------------------------------------------------------------------------------- %    
    gamma = 1; beta = 1; c = model.const.c;
    acc_genRel = model.sat.rel.*model.centralPlanet.GM/(c^2*norm(xJ2000)^3)*((2*(beta+gamma)*model.centralPlanet.GM/norm(xJ2000)-...
        gamma*norm(vJ2000)^2).*xJ2000+2*(1+gamma)*xJ2000'*vJ2000.*vJ2000);
    fprintf(TextAcc, "RelatCorr   "+num2str(norm(acc_genRel))+"\n");
    % ---------------------------------------------------------------------------------------------------------- %

    % acceleration due to the solar radiation pressure of Sun
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_SRPSun = accelsrp(xJ2000,model.sat.srp,model.const,'SUN',t,model.frame.integr,model.centralPlanet.stringName,model);
    fprintf(TextAcc, "SunRadPress   "+num2str(norm(acc_SRPSun))+"\n");
    % ---------------------------------------------------------------------------------------------------------- %
    
    % acceleration due to the Earth albedo
    % ---------------------------------------------------------------------------------------------------------- %    
    acc_alb = accelalb(xJ2000,model.sat,model.const,'EARTH',t,model.frame.integr,model.centralPlanet.stringName,model);
    fprintf(TextAcc, "EarthAlb   "+num2str(norm(acc_alb))+"\n");
    % ---------------------------------------------------------------------------------------------------------- %

    % equations of motion
    % ---------------------------------------------------------------------------------------------------------- %
    xdot   = vJ2000;
    vdot   = acc_centralPlanet + acc_earth + acc_pointMasses + acc_SRPSun + acc_alb + acc_genRel;
    % ---------------------------------------------------------------------------------------------------------- %

    % output
    % ---------------------------------------------------------------------------------------------------------- %
    X = [xdot;vdot];
    % ---------------------------------------------------------------------------------------------------------- %   
    
    if isfield(model,'wb')
        wb = model.wb;
        waitbar((t-wb.T(1))/(wb.T(end)-wb.T(1)),wb.bar,sprintf('propagating... %07.3f %%',100*(t-wb.T(1))/(wb.T(end)-wb.T(1))));
    end
end
