function [x,y,lat,lon,thk,topg,usrf] = buildant(doplot,filename)
% BUILDANT  helper function for ant.m, to build a Matlab/Octave set of
% variables, for Antarctic ice sheet, on 50km grid, from existing NetCDF file
% example which plots 3 fields:
%   >> buildant
% example which fills these variables, but no plot
%   >> [x,y,lat,lon,thk,topg,usrf] = buildant(0);
% as above but attempting to read a different NetCDF file:
%   >> [x,y,lat,lon,thk,topg,usrf] = buildant(0,'foobar.nc');

if nargin < 1, doplot = 1; end
if nargin < 2, filename = 'Antarctica_5km_dev1.0.nc'; end

disp(['reading variables x,y,lat,lon,prcp,thk,topg,usrf from NetCDF file ' filename])

S = netcdf(filename);  % reads NetCDF file into a large structure

x = double(S.VarArray(14).Data);
y = double(S.VarArray(15).Data);
lat = squeeze(double(S.VarArray(2).Data));
lon = squeeze(double(S.VarArray(3).Data));
thk = squeeze(double(S.VarArray(25).Data));
topg = squeeze(double(S.VarArray(21).Data));
usrf = squeeze(double(S.VarArray(22).Data));

if doplot==0, return; end

disp('plotting 2 fields')
figure(1)
surf(x/1000,y/1000,usrf), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('surface elevation "usrf"  (m)'), colorbar
figure(2)
surf(x/1000,y/1000,topg), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('bed elevation "topg"  (m)'), colorbar
figure(3)

