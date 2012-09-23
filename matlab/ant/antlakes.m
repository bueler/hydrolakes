function antlakes(tyears,nn)
% ANTLAKES  Use Antarctic geometry data from SeaRISE, and assumed uniform
% melt rate, to compute water distribution solving a water conservation
% model which assumes that pressure is overburden.
% form:  antlakes(tyears,nn)
% where    tyears = run time in years [default=100.0]
%              nn = stride for subsample [default=10 for 50 km res]
% example:  >> antlakes(100.0,10)
% Calls:  CONSERVE, BUILDANT
% Depends: NETCDF

if nargin<1, tyears=100.0; end
if nargin<2, nn=10; end

filename = 'Antarctica_5km_dev1.0.nc';
fprintf('reading variables x,y,lat,lon,thk,topg,usrf\n  from NetCDF file %s\n',filename)
[x,y,lat,lon,thk,topg,usrf] = buildant(0,filename);

fprintf('x range: %.3f km to %.3fkm\n',min(x)/1000,max(x)/1000)
fprintf('y range: %.3f km to %.3fkm\n',min(y)/1000,max(y)/1000)

dx = x(2)-x(1);
fprintf('subsampling to %d km resolution\n', dx*nn / 1000.0)
clear dx
x = x(1:nn:end);
y = y(1:nn:end);
lat = lat(1:nn:end,1:nn:end);
lon = lon(1:nn:end,1:nn:end);
thk = thk(1:nn:end,1:nn:end);
topg = topg(1:nn:end,1:nn:end);
usrf = usrf(1:nn:end,1:nn:end);

fprintf('showing topg, usrf, floatmask\n')

figure(1)
imagesc(x/1000,flipud(-y)/1000,flipud(topg)), colorbar
title('bed elevation')
xlabel('x (km)'), ylabel('y (km)')

figure(2)
imagesc(x/1000,flipud(-y)/1000,flipud(usrf)), colorbar
title('surface elevation')
xlabel('x (km)'), ylabel('y (km)')

floatmask = (usrf - topg > thk+10.0);
figure(3)
imagesc(x/1000,flipud(-y)/1000,flipud(floatmask)), colorbar
title('mask for where ice is floating')
xlabel('x (km)'), ylabel('y (km)')

spera = 31556926.0;

Phi0 = 0.01 / spera;   % 1 cm a-1
Phi = zeros(size(floatmask));
Phi(floatmask < 0.5) = Phi0;

W0 = zeros(size(thk));
ts = 0.0;
te = tyears*spera;
W = conserve(x,y,topg,usrf,floatmask,W0,Phi,ts,te,5);

figure(4)
imagesc(x/1000,flipud(-y)/1000,flipud(W)), colorbar
title(sprintf('water amount from CONSERVE after time of %.3f year',(te-ts)/spera))
xlabel('x (km)'), ylabel('y (km)')

