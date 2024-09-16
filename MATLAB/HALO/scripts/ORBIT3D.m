% ORBIT3D allows the view of the trajectory of a lunar orbiter in a 3D window.
% please make sure the mice toolbox is installed, otherwise download it at the following link:
% https://naif.jpl.nasa.gov/naif/toolkit_MATLAB_PC_Windows_VisualC_MATLAB8.x_64bit.html
%
% In the 3D window the user can visualize the motion of the orbiter, along with the ground track,
% the orbit track, the Sun direction (yellow vector), the Earth direction (green vector), and all
% the reference frames. In particular:
% ECI       --> J2000 frame
% SCSF      --> MOON_ME frame
% X0,Y0,Z0  --> SCSF at t = t0
% For further details on the frames please see the NAIF definition.
%
% The instructions included in the lines %%%%%% can be edited.
% - ECI_frame   = 1 to view the J2000 frame, 0 otherwise
% - SCSF_frame  = 1 to view the Surface Centered Surface Fized frame, 0 otherwise
% - gtk_on      = 1 to view the ground track in the 3D window, 0 otherwise
% - planisphere = 1 to view the 2D window with the ground track, 0 otherwise
% - orbittrack  = 1 to view the orbit track in the 3D window, 0 otherwise
% - playingSpeed = [double] to play the animation at different velocity (10 suggested)
% - fact         = [double] the half dimension of the equatorial plane (2000 suggested)
%
% Note: work is in progress to improve the code. 
%
% (c)Copyright 2017 - Ennio Condoleo
% Sapienza University of Rome, 08-25-2017
% ennio.condoleo@uniroma1.it

%% INIT. WORK ENVIRONMENT
% =========================================================================================================================================================
    % clear workspace and kernel pool
    clearvars; clear; clc;
    cspice_kclear;
    % current directory
    pathRun = cd;
    % adding tools: lunarhpop.m, dragzoom.m, force_delete.m
    pathToTools = [pathRun,'/tools'];
    addpath(pathToTools);
    % loading spice kernels
    cspice_furnsh('metakernel.tm');
    % loading ORB data
    load('output/ORBdata.mat');
% =========================================================================================================================================================

%% USER'S INPUT
% =========================================================================================================================================================
    % flags
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ECI_frame   = 1;   % [1/0] turn on to view the ECI  frame
    SCSF_frame  = 1;   % [1/0] turn on to view the SCSF frame
    gtk_on      = 1;   % [1/0] turn on to view the ground track
    planisphere = 1;   % [1/0] turn on to view the planisphere
    orbittrack  = 1;   % [1/0] turn on to view the orbit track in the 3D window
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % setting for the video
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    playingSpeed = 10;    % speed of the animation (suggested 10)
    fact         = 3000;  % equatorial plane [km] (it should be greater than 1738 km)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================================================================================================

%% INIT. VARIABLES
% =========================================================================================================================================================
    % creating ephemerides file for the orbiter
    createBSPfile(orb);
    cspice_furnsh('metakernel.tm');

    lengthaxis  = norm(orb.XJ2000(1,1:3));
    image_file  = 'input/Moon_grid.bmp';
    RE          = 1737.74;
    R_ECI_SCSF  = cspice_pxform('J2000','MOON_ME',orb.t');
    R_ECI_INER  = cspice_pxform('J2000','MOON_ME',orb.t(1));
    R_SCSF_INER = cspice_pxfrm2('MOON_ME','MOON_ME',orb.t',orb.t(1).*ones(size(orb.t')));
    count       = 0;
    step        = max([floor(playingSpeed/(orb.t(2)-orb.t(1))),1]);

    % Sun and Earth direction
    t     = orb.epoch.span;
    Xs = cspice_spkezr('SUN',t,'J2000','NONE','MOON');
    Xe = cspice_spkezr('EARTH',t,'J2000','NONE','MOON');
    for j = 1:size(Xs,2)
        Xs(1:3,j) = Xs(1:3,j).*fact./norm(Xs(1:3,j));
        Xe(1:3,j) = Xe(1:3,j).*fact./norm(Xe(1:3,j));
        Xs(1:3,j) = R_ECI_INER*Xs(1:3,j);
        Xe(1:3,j) = R_ECI_INER*Xe(1:3,j);
    end
% =========================================================================================================================================================

%% FIGURE
% =========================================================================================================================================================    
    % 3D animation
    f1 = figure('Name','Orbit3D');
    set(f1,'Color','k','Units','normalized',...
        'Position',[0.507291666666667 0.433333333333333 0.490104166666667 0.49]);
    set(f1,'Renderer','opengl');
    ax1 = axes;
    ax1.Color = [0,0,0];
    ax1.Visible = 'off'; % axis off
    hold on; material dull; lighting gouraud; l = light; 
    ax1.DataAspectRatio = [1,1,1]; % axis equal;

    
    % Inertial frame
    plot3(ax1,0,0,0,'w*','MarkerSize',10); grid on;
    iner(1) = plot3(ax1,[0,lengthaxis],[0,0],[0,0],'w-','LineWidth',2);
    iner(2) = plot3(ax1,[0,0],[0,lengthaxis],[0,0],'w-','LineWidth',2);
    iner(3) = plot3(ax1,[0,0],[0,0],[0,lengthaxis],'w-','LineWidth',2);
    iner(4) = plot3(ax1,[0,0],[0,0],[0,-lengthaxis],'w--','LineWidth',2);
    iner(5) = plot3(ax1,[0,0],[0,-lengthaxis],[0,0],'w--','LineWidth',2);
    iner(6) = plot3(ax1,[0,-lengthaxis],[0,0],[0,0],'w--','LineWidth',2);
    iner(7) = text(lengthaxis,0,0,'X_0','Color','w');
    iner(8) = text(0,lengthaxis,0,'Y_0','Color','w');
    iner(9) = text(0,0,lengthaxis,'Z_0','Color','w');
    
    % Moon
    cdata = imread(image_file);view(175,60);
    [x,y,z] = ellipsoid(ax1,0, 0, 0, RE, RE, RE, 100);
    globe   = surface(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    hg = hgtransform('Parent',ax1);
    set(globe,'Parent',hg,'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    
    % SCSF frame
    if SCSF_frame
        axis_SCSF(1) = plot3(ax1,[0,lengthaxis],[0,0],[0,0],'b','LineWidth',2);
        axis_SCSF(2) = plot3(ax1,[0,0],[0,lengthaxis],[0,0],'b','LineWidth',2);
        axis_SCSF(3) = plot3(ax1,[0,0],[0,0],[0,lengthaxis],'b','LineWidth',2);
        tx_SCSF(1)   = text(lengthaxis,0,0,'X_S_C_S_F','Color','b');
        tx_SCSF(2)   = text(0,lengthaxis,0,'Y_S_C_S_F','Color','b');
        tx_SCSF(3)   = text(0,0,lengthaxis,'Z_S_C_S_F','Color','b');
        set(axis_SCSF,'Parent',hg);
        set(tx_SCSF,'Parent',hg);
    end
    
    % ECI frame
    if ECI_frame
        axis_ECI(1) = plot3(ax1,[0,lengthaxis*R_ECI_INER(1,1)],[0,lengthaxis*R_ECI_INER(2,1)],[0,lengthaxis*R_ECI_INER(3,1)],'c','LineWidth',2);
        axis_ECI(1) = plot3(ax1,[0,lengthaxis*R_ECI_INER(1,2)],[0,lengthaxis*R_ECI_INER(2,2)],[0,lengthaxis*R_ECI_INER(3,2)],'c','LineWidth',2);
        axis_ECI(1) = plot3(ax1,[0,lengthaxis*R_ECI_INER(1,3)],[0,lengthaxis*R_ECI_INER(2,3)],[0,lengthaxis*R_ECI_INER(3,3)],'c','LineWidth',2);
        tx_ECI(1)   = text(lengthaxis*R_ECI_INER(1,1),lengthaxis*R_ECI_INER(2,1),lengthaxis*R_ECI_INER(3,1),'X_E_C_I','Color','c');
        tx_ECI(1)   = text(lengthaxis*R_ECI_INER(1,2),lengthaxis*R_ECI_INER(2,2),lengthaxis*R_ECI_INER(3,2),'Y_E_C_I','Color','c');
        tx_ECI(1)   = text(lengthaxis*R_ECI_INER(1,3),lengthaxis*R_ECI_INER(2,3),lengthaxis*R_ECI_INER(3,3),'Z_E_C_I','Color','c');
    end
    
    % equatorial plane
    [xp,yp] = meshgrid(-fact:fact:fact);
    equat = surf(ax1,xp,yp,zeros(size(xp)));
    equat.FaceColor = [0.9,0.9,0.9];
    equat.EdgeColor = [0.9,0.9,0.9];
    equat.LineWidth = 0.1;
    equat.FaceLighting = 'none';
    equat.EdgeLighting = 'none';
    equat.FaceAlpha = 0.5;
    equat.EdgeAlpha = 0.5;
    
    % enable muose commands to control the figure
    % see help dragzoom
    try
        dragzoom(f1,'on');
        str = {     'MOUSE CONTROLS';...
                    'hold left-button: pan';...
                    'hold right-button: rotate';...
                    'scroll mean-button: zoom'};
        % annotation mouse        
        ann.mouse = annotation(f1,'textbox',[0.8,0.88,0.1,0.1],'String',str,'Color','w',...
            'FitBoxToText','on','EdgeColor',[1,1,1],'LineWidth',1.5,...
            'FontWeight','normal','FontSize',11);
    catch
    end

    % annotation time
    str2 = {     'TIME';...
                 'SMA       [km]: ';...
                 'ECC           : ';...
                 'INC      [deg]: ';...
                 'LAN      [deg]: ';...
                 'AOP      [deg]: ';...
                 'MA       [deg]: '};  
    ann.time = annotation(f1,'textbox',[0.01,0.86,0.1,0.1],'String',str2,'Color','w',...
            'FitBoxToText','on','EdgeColor',[1,1,1],'LineWidth',1.5,...
            'FontWeight','normal','FontSize',11);
    
    % rewind button
    rew=uicontrol('Parent',f1,'Style','pushbutton','String','Rewind',...
    'Units','normalized','Position',[0.89 0.6 0.1 0.08],'Visible','on','Callback',...
    'set(rew, ''UserData'', rand);repl=1;');
    rew.BackgroundColor = [0.1,0.1,0.1];
    rew.ForegroundColor=[1,1,1];
    rew.FontSize = 12;
    rew.FontWeight = 'bold';
    rew.DeleteFcn = 'rmpath(pathToTools)';
    
    % computing trajectory
    cspice_boddef('orbiter',-1000001);
    cspice_furnsh('orbEph.bsp');
    indeces = 1:step:length(t);
    Xo_ECI = cspice_spkezr('orbiter',t,'J2000','NONE','MOON');
    X      = zeros(size(Xo_ECI));
    X_SCSF = zeros(size(Xo_ECI));
    for j = 1:size(Xo_ECI,2)
        X(1:3,j) = R_ECI_INER*Xo_ECI(1:3,j);
        X(4:6,j) = R_ECI_INER*Xo_ECI(4:6,j);
        X_SCSF(1:3,j) = R_ECI_SCSF(:,:,j)*Xo_ECI(1:3,j);
    end
    if orbittrack
        plot3(X(1,:),X(2,:),X(3,:),'r--','LineWidth',1);
    end
    % keplerian elements
    elts = cart2kepl(X,t);
    % latitude and longitude
    [lat_,lon_,alt_] = scsf2lla(X(1:3,:),RE,0.0012);
    [lat,lon,alt] = scsf2lla(X_SCSF(1:3,:),RE,0.0012);
    if gtk_on
        gtk = plot3(ax1,RE.*cos(lat).*cos(lon),...
                  RE.*cos(lat).*sin(lon),...
                  RE.*sin(lat),'Color',[190,74,34]./255,'LineStyle','-',...
                  'MarkerFaceColor',[190,74,34]./255,'MarkerSize',1,'LineWidth',1);
        set(gtk,'Parent',hg);
    end
    
    % planisphere
    if planisphere
        f2 = figure('Name','planisphere');
        set(f2,'Color','w','Units','normalized',...
        'Position',[0.007291666666667 0.433333333333333 0.490104166666667 0.49]);
        set(f2,'Renderer','opengl');
        
        hold on;
        LongSSP = lon.*180/pi; 
        LatSSP  = lat.*180/pi;  
        LongSSP(LongSSP>180) = LongSSP(LongSSP>180)-360;
        set(gca,'XTick',[-180 -150 -120 -90 -60 -30 0 ...
            30 60 90 120 150 180],'XTickMode','manual');
        set(gca,'YTick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'YTickMode','manual');
        imagesc([-180,180],[90,-90],cdata);
        grid on;
        axis([-180 180 -90 90]);
        plot(LongSSP,LatSSP,'Color',[190,74,34]./255,'Marker','.','MarkerSize',5,...
            'LineStyle','none','MarkerFaceColor',[190,74,34]./255,'MarkerEdgeColor',[190,74,34]./255);
        ax2 = gca;
    end
    
    % animation starts
    while 1
        try
            for j = indeces       
                count = count+1;
                set(hg,'Matrix',[R_SCSF_INER(:,:,j),zeros(3,1);zeros(1,3),1]);
                p(1) = plot3(ax1,[0,Xs(1,j)],[0,Xs(2,j)],[0,Xs(3,j)],'y-','LineWidth',1.5);
                p(2) = plot3(ax1,[0,Xe(1,j)],[0,Xe(2,j)],[0,Xe(3,j)],'g-','LineWidth',1.5);
                p(3) = plot3(ax1,X(1,j),X(2,j),X(3,j),'ro','MarkerFaceColor','r','MarkerSize',5);
                if gtk_on
                    p(4) = plot3(ax1,RE*cos(lat_(j))*cos(lon_(j)),...
                        RE*cos(lat_(j))*sin(lon_(j)),...
                        RE*sin(lat_(j)),'Color',[190,74,34]./255,'Marker','square',...
                        'MarkerFaceColor',[190,74,34]./255,'MarkerEdgeColor','w','MarkerSize',5);
                end
                % planisphere animation
                if planisphere
                    p(5) = plot(ax2,LongSSP(j),LatSSP(j),'Color',[190,74,34]./255,'Marker','square','MarkerSize',8,...
            'MarkerFaceColor',[190,74,34]./255,'MarkerEdgeColor','w');
                end
                % updating annotation time 
                epoch_time = cspice_et2utc(t(j),'C',3);
                ann.time.String = { sprintf('TIME: %s',epoch_time);...
                                    sprintf('SMA [km]: %10.2f',elts(1,j));...
                                    sprintf('ECC: %8.5f',elts(2,j));...
                                    sprintf('INC [deg]: %7.3f',elts(3,j));...
                                    sprintf('LAN [deg]: %7.3f',elts(4,j));...
                                    sprintf('AOP [deg]: %7.3f',elts(5,j));...
                                    sprintf('MA [deg]: %7.3f',elts(6,j))}; 

                l.Position = Xs(1:3,j);
                drawnow;
                if j<length(t)
                    delete(p);
                end
            end
        catch
            % cleaning work environment
            cspice_kclear;
            fclose('all');
            close all;
            delete('orbEph.bsp');
            rmpath(pathToTools);
            clc;
            return;
        end
        waitfor(rew,'UserData');            
        count = 0;
        delete(p);
    end
% =========================================================================================================================================================

