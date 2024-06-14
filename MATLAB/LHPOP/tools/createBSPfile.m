
%% Creating a SPK ephemeris file to export the results (You can use it in STK free version)
% Define the segment identifier parameters
function [] = createBSPfile(orb)
try   
   body       = -1000001;
   cspice_boddef('lunar_orbiter',body);
   center     = cspice_bodn2c(orb.centralPlanet.stringName);
   ref        = 'J2000';
   poly_deg   = 9;
   spk8       = 'orbEph.bsp';
   segid = 'SPK type 8 test segment';         % Create a segment identifier.
   if (exist(spk8,'file')==2)
      force_delete(spk8);        % Open a new SPK file. Delete if a file of the same name exists.
   end
   spkhandle = cspice_spkopn(spk8,segid,4); 
   cspice_spkw08(spkhandle,body,center,ref,orb.t(1),orb.t(end),segid,poly_deg,orb.XJ2000',orb.t(1),orb.epoch.span(2)-orb.epoch.span(1));        % Create a type 8 segment.
   cspice_spkcls(spkhandle);        % Close the SPK file.
catch
   spk8       = regexprep([date,'_',num2str(rem(floor(abs(rand*1e4)),9999)),'.bsp'],'-','_');
   spkhandle = cspice_spkopn(spk8,segid,4); 
   cspice_spkw08(spkhandle,body,center,ref,orb.t(1),orb.t(end),segid,poly_deg,orb.XJ2000',orb.t(1),orb.epoch.span(2)-orb.epoch.span(1));        % Create a type 8 segment.
   cspice_spkcls(spkhandle);        % Close the SPK file.
    fprintf([' \n bsp file with a random name generated!',...
        '\n This error occurs when the permission to delete spk8.bsp is denied',...
        '\n The random named file will have to be deleted (if possible) at the end of the run\n',...
        'to avoid the storage of many files\n']);
    fprintf(['file generated: %s\n',spk8]);
end
cspice_kclear;
end
