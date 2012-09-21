function [x,y,thk,topg,usurf,icemask, outline] = buildnbreen(doplot,filename)
% BUILDNBREEN  helper function for NBREENWATER, to build a Matlab/Octave set of
% variables from NetCDF file nbreen_input.nc
% example which plots 3 fields:
%   >> buildant
% example which fills these variables, but no plot
%   >> [x,y,lat,lon,thk,topg,usrf] = buildant(0);
% as above but attempting to read a different NetCDF file:
%   >> [x,y,lat,lon,thk,topg,usrf] = buildant(0,'foobar.nc');

if nargin < 1, doplot = 1; end
if nargin < 2, filename = 'nbreen_input.nc'; end

S = netcdf(filename);  % reads NetCDF file into a large structure

y = double(S.VarArray(1).Data);
x = double(S.VarArray(2).Data);
topg = squeeze(double(S.VarArray(3).Data));
thk = squeeze(double(S.VarArray(4).Data));
usurf = squeeze(double(S.VarArray(5).Data));
icemask = squeeze(double(S.VarArray(6).Data));
outline = squeeze(double(S.VarArray(7).Data));
for i=1:size(icemask,1)
    for j=1:size(icemask,2)
%         if (i==1 || j==1 || i==size(icemask,1) || j==size(icemask,2))
%             icemask(i,j) = 0;
%         end
        if usurf(i,j)==0.0, outline(i,j)=0; end
    end
end

if doplot==0, return; end

disp('plotting 2 fields')
figure(1)
surf(y/1000,x/1000,usurf), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('surface elevation "usurf"  (m)'), colorbar
figure(2)
surf(y/1000,x/1000,topg), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('bed elevation "topg"  (m)'), colorbar

