function [W, P] = nbreenwater(tyears,nn,W0,P0)
% NBREENWATER  **FIXME this comment**
% Form:  [W, P] = nbreenwater(tyears,nn,W0,P0)
% All input and output arguments are optional (and have given defaults).
% Calls:   BUILDNBREEN, CONSERVEWATER, DOUBLEDIFF
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

p = params();

Phi0arr = zeros(size(icemask));
Phi = zeros(size(icemask));
for i=1:length(x)
    for j=1:length(y)
        Phi0arr(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/p.spera; 
        if outline(i,j)>0.5, Phi(i,j) = Phi0arr(i,j); end
    end
end

if true
  fprintf('showing initial fields: topg, usurf, outline, Phi, geometric psi\n')

  figure(1)
  set(gcf,'position',[100 400 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(topg)), colorbar
  title('bed elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(usurf)), colorbar
  title('surface elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(2)
  set(gcf,'position',[100 300 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(outline)), colorbar
  title('outline')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(Phi*p.spera)), colorbar
  title('water input  (m/a)')
  xlabel('x (km)'), ylabel('y (km)')

end

if nargin < 3
  W0 = zeros(size(topg));
end
if nargin < 4
  P0 = zeros(size(topg));
end

vb = 50.0/p.spera;    % ice speed (m s-1; = 50 m a-1)
magvb = vb * ones(size(topg));
ts = 0.0;
te = tyears*p.spera;
fprintf('calling doublediff() to do run for %.3f years:\n\n',tyears)
%W = conservewater(x,y,topg,usurf,outline,W0,Phi,ts,te,5);
%[W, Y, P] = damper(x,y,topg,usurf,magvb,outline,W0,Y0,Phi,ts,te,5);
[W, P] = doublediff(x,y,topg,usurf,magvb,outline,W0,P0,Phi,ts,te,5);

if true
  fprintf('showing final fields: W, P\n')

  figure(3)
  set(gcf,'position',[100 200 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(W)), colorbar
  title('water thickness W  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  Pmask = ones(size(P));
  Po = p.rhoi * p.g * thk;
  Pmask(P < 0.001*Po) = 0;
  Pmask(P > 0.999*Po) = 2;
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(Pmask)), colorbar
  title('pressure mask:  0=under, 1=normal, 2=over')
  xlabel('x (km)'), ylabel('y (km)')

  figure(4)
  set(gcf,'position',[100 100 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(Po),[0 5.0e6]), colorbar
  title('overburden pressure P_o  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(P),[0 5.0e6]), colorbar
  title('water pressure P  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(5)
  Wpos=W(W>0);  Ppos=P(W>0);  Popos=Po(W>0);
  plot(Wpos(:),Ppos(:)./Popos(:),'ko','markersize',8)
  hold on, WFC=0:.01:1.5*p.Wr; plot(WFC,(WFC/p.Wr).^(7/2),'r--'), hold off
  axis([0 1.5*p.Wr 0 1.2])
  xlabel('water thickness W  (m)')
  ylabel('water pressure P as fraction of overburden')
  title('is P=P(W)?  (Flowers & Clarke is red dashed)')
end

