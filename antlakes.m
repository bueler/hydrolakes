% ANTLAKES  FIXME:  for now it just plots fields

filename = 'Antarctica_5km_dev1.0.nc';
fprintf('reading variables x,y,lat,lon,thk,topg,usrf from NetCDF file %s\n',filename)
[x,y,lat,lon,thk,topg,usrf] = buildant(0,filename);

fprintf('subsampling to 25km resolution\n')
x = x(1:5:end);
y = y(1:5:end);
lat = lat(1:5:end,1:5:end);
lon = lon(1:5:end,1:5:end);
thk = thk(1:5:end,1:5:end);
topg = topg(1:5:end,1:5:end);
usrf = usrf(1:5:end,1:5:end);

fprintf('showing topg, usrf, floatmask\n')

figure(1)
imagesc(x/1000,flipud(-y)/1000,flipud(topg)), view(2), colorbar
title('bed elevation')
xlabel('x (km)'), ylabel('y (km)')

figure(2)
imagesc(x/1000,flipud(-y)/1000,flipud(usrf)), view(2), colorbar
title('surface elevation')
xlabel('x (km)'), ylabel('y (km)')

floatmask = (usrf - topg > thk+10.0);
figure(3)
imagesc(x/1000,flipud(-y)/1000,flipud(floatmask)), view(2), colorbar
title('mask for where ice is floating')
xlabel('x (km)'), ylabel('y (km)')

spera = 31556926.0;
W0 = zeros(size(thk));
ts = 0.0;
te = spera;
Wnew = conserve(x,y,topg,usrf,floatmask,W0,ts,te,10);
figure(4)
imagesc(x/1000,flipud(-y)/1000,flipud(Wnew)), view(2), colorbar
title(sprintf('water amount from CONSERVE after time of %f sec',te-ts))
xlabel('x (km)'), ylabel('y (km)')

