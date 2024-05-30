function [lat,lon,alt] = scsf2lla(varargin)

if(length(varargin)==3)
    x=varargin{1}(1,:);
    y=varargin{1}(2,:);
    z=varargin{1}(3,:);
    req_Surface=varargin{2};
    e_Surface=varargin{3};
elseif(length(varargin)==5)
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    req_Surface=varargin{4};
    e_Surface=varargin{5};
else
    disp('Error: usage [lat,lon,alt] = scsf2lla(x,y,z,req_Surface,e_Surface)');
    disp('or                         = scsf2lla(r_SCSF,req_Surface,e_Surface)');
    return;
end

a = req_Surface;
e2 = e_Surface^2;

% calculations:
b   = sqrt(a^2*(1-e2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);
lon = atan2(y,x);
lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e2.*a.*cos(th).^3));
N   = a./sqrt(1-e2.*sin(lat).^2);
alt = p./cos(lat)-N;

% return lon in range [0,2*pi)
lon = mod(lon,2*pi);

% correct for numerical instability in altitude near exact poles:
% (after this correction, error is about 2 millimeters, which is about
% the same as the numerical precision of the overall function)

k=abs(x)<1 & abs(y)<1;
alt(k) = abs(z(k))-b;

return
