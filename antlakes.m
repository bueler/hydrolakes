function antlakes(tyears,nn)
% ANTLAKES  FIXME:  for now it just plots fields

if nargin<1, tyears=100.0; end
if nargin<2, nn=10; end

filename = 'Antarctica_5km_dev1.0.nc';
fprintf('reading variables x,y,lat,lon,thk,topg,usrf from NetCDF file %s\n',filename)
[x,y,lat,lon,thk,topg,usrf] = buildant(0,filename);

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
W0 = zeros(size(thk));
ts = 0.0;
te = tyears*spera;
[W,xx,yy,qx,qy] = conserve(x,y,topg,usrf,floatmask,W0,ts,te,10);

figure(4)
imagesc(x/1000,flipud(-y)/1000,flipud(W)), colorbar
title(sprintf('water amount from CONSERVE after time of %.3f year',(te-ts)/spera))
xlabel('x (km)'), ylabel('y (km)')

return

figure(5), scale=1.0e11;
quiver(xx/1000,yy/1000,scale*qx,scale*qy)
title('ice flux  (m^2 s^-1)'), xlabel('x (km)'), ylabel('y (km)')

