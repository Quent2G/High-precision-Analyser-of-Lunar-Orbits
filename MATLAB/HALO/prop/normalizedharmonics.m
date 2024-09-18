% function [Cnm,Snm] = normalizedharmonics(filepath,maxdegree)
% 
% DESCRIPTION
% --------------------------------------------------------------------------
% normalizedharmonics read the normalized coefficient of the harmonics
% from filepath. The file of the harmonics' coefficients has to be organised as:
% 0   0    1.0 0.0
% 1   0    0.0 0.0
% 1   1    0.0 0.0
% 2	0	-9.0933087157900000E-5	0.0000000000000000E+00
% 2	1	-2.7220323615900000E-09	-7.5751829208300000E-10
% 2	2	3.4709851601400000E-5	1.6729490538300000E-08
% 3	0	-3.2035914003000000E-06	0.0000000000000000E+00
% 3	1	2.6327440121800000E-05	5.4643630898200000E-06
% 3	2	1.4188179329400000E-05	4.8920365004800000E-06
% ...
% where the first two columns are related to the degree and order of the harmonics'
% coefficients.
%
% INPUT 
% --------------------------------------------------------------------------
% filepath  = path to file of the harmonics' coefficients
% maxdegree = maximum degree for reading the file of harmonics
%
% OUTPUT
% --------------------------------------------------------------------------
% Cnm      = coefficient of the cosine-term (nm) in the gravity potential
% Snm      = coefficient of the sine-term (nm) in the gravity potential
%
% AUTHOR
% --------------------------------------------------------------------------
% Ennio Condoleo,
% Jan 02, 2017 - Rome
% ennio.condoleo@uniroma1.it
%
% See also accelharmonic prophpop
%
function [Cnm,Snm] = normalizedharmonics(filepath,maxdegree)

    Cnm = zeros(maxdegree+1,maxdegree+1);
    Snm = zeros(maxdegree+1,maxdegree+1);
    fid = fopen(filepath,'r');
    for n=0:maxdegree
        for m=0:n
            temp = fscanf(fid,'%d %d %f %f',[4 1]);        
            Cnm(n+1,m+1) = temp(3);
            Snm(n+1,m+1) = temp(4);
        end
    end
    fclose(fid);
end
