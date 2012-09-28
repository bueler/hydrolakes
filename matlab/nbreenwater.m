function [W, Y] = nbreenwater(tyears,nn,W0,Y0)
% NBREENWATER  **FIXME this comment**
% Form:  W = nbreenwater(tyears,nn)
% (Note all input and output arguments are optional.)
% Calls:   BUILDNBREEN, CONSERVEWATER, DAMPER
% Depends: NETCDF

if nargin<1, tyears=5.0; end
if nargin<2, nn=4; end

filename = 'nbreen_input.nc';

fprintf('reading variables x,y,thk,topg,usurf,icemask\n  from text files %s\n',filename)
[x,y,thk,topg,usurf,icemask,outline] = buildnbreen(0,filename);

fprintf('x range: %.3f km to %.3fkm\n',min(x)/1000,max(x)/1000)
fprintf('y range: %.3f km to %.3fkm\n',min(y)/1000,max(y)/1000)

dx = x(2)-x(1);
fprintf('subsampling to %d m resolution\n', dx*nn)
clear dx
x = x(1:nn:end);
y = y(1:nn:end);
thk = thk(1:nn:end,1:nn:end);
topg = topg(1:nn:end,1:nn:end);
usurf = usurf(1:nn:end,1:nn:end);
icemask = icemask(1:nn:end,1:nn:end);
outline = outline(1:nn:end,1:nn:end);

spera = 31556926.0;
Phi0arr = zeros(size(icemask));
Phi = zeros(size(icemask));
for i=1:length(x)
    for j=1:length(y)
        Phi0arr(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/spera; 
        if outline(i,j)>0.5, Phi(i,j) = Phi0arr(i,j); end
        if usurf(i,j)==0.0, usurf(i,j)=topg(i,j); end
    end
end

if true
%if nargout < 1
  fprintf('showing initial fields: topg, usurf, outline, Phi, geometric psi\n')

  figure(1)
  imagesc(x/1000,flipud(-y)/1000,flipud(topg)), colorbar
  title('bed elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  %set(gcf,'PaperPositionMode','auto'), print('-dpsc2',strcat('fig_bed.eps'))

  figure(2)
  imagesc(x/1000,flipud(-y)/1000,flipud(usurf)), colorbar
  title('surface elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  %set(gcf,'PaperPositionMode','auto'), print('-dpsc2',strcat('fig_usurf.eps'))

  figure(3)
  imagesc(x/1000,flipud(-y)/1000,flipud(outline))
  title('outline')
  xlabel('x (km)'), ylabel('y (km)')
  %set(gcf,'PaperPositionMode','auto'), print('-dpsc2',strcat('fig_outline.eps'))

  figure(4)
  imagesc(x/1000,flipud(-y)/1000,flipud(Phi*spera)), colorbar
  title('water input  (m/a)')
  xlabel('x (km)'), ylabel('y (km)')
  %set(gcf,'PaperPositionMode','auto'), print('-dpsc2',strcat('fig_melt.eps'))

  figure(5)
  rhoi = 910.0;  rhow = 1028.0;  g = 9.81;
  psi = rhoi * g * thk + rhow * g * topg;
  imagesc(x/1000,flipud(-y)/1000,flipud(psi)), colorbar
  title('hydraulic potential from geometry  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')
end

%return

if nargin < 3
  W0 = zeros(size(topg));
end
if nargin < 4
  Y0 = W0;
end

vb = 50.0/spera;    % ice speed (m s-1; = 50 m a-1)
magvb = vb * ones(size(topg));
ts = 0.0;
te = tyears*spera;
%W = conservewater(x,y,topg,usurf,outline,W0,Phi,ts,te,5);
[W, Y, P] = damper(x,y,topg,usurf,magvb,outline,W0,Y0,Phi,ts,te,5);

if true
%if nargout < 1
  fprintf('showing final fields: W, Y, P\n')

  figure(6)
  imagesc(x/1000,flipud(-y)/1000,flipud(W)), colorbar
  title('water thickness W  (m)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(7)
  imagesc(x/1000,flipud(-y)/1000,flipud(Y)), colorbar
  title('capacity thickness Y  (m)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(8)
  imagesc(x/1000,flipud(-y)/1000,flipud(P)), colorbar
  title('water pressure P  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(9)
  Pmask = ones(size(P));
  Po = rhoi * g * thk;
  Pmask(P < 0.001*Po) = 0;
  Pmask(P > 0.999*Po) = 2;
  imagesc(x/1000,flipud(-y)/1000,flipud(Pmask)), colorbar
  title('pressure mask:  0 = underpressure,  1 = normal,  2 = overpressure')
  xlabel('x (km)'), ylabel('y (km)')
end

