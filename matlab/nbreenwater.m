function [x, y, W, P] = nbreenwater(tyears,nn,W0,P0)
% NBREENWATER  **FIXME this comment**
% Form:  [x, y, W, P] = nbreenwater(tyears,nn,W0,P0)
% All input and output arguments are optional (and have given defaults).
% Calls:   BUILDNBREEN, CONSERVEWATER, DOUBLEDIFF
% Depends: NETCDF

if nargin<1, tyears=5.0; end
if nargin<2, nn=4; end

filename = 'nbreen_input.nc';

fprintf('reading variables x,y,thk,topg,usurf,icemask\n  from text files %s\n',filename)
[x,y,~,topg,usurf,~,outline] = buildnbreen(0,filename);

dx = x(2)-x(1);
fprintf('subsampling to %d m resolution\n', dx*nn)
x = x(1:nn:end);
y = y(1:nn:end);
topg = topg(1:nn:end,1:nn:end);
usurf = usurf(1:nn:end,1:nn:end);
outline = outline(1:nn:end,1:nn:end);

% consistent thickness (ignores values in file)
thk = usurf - topg;
thk(thk < 0) = 0.0;

% create floating and ice free masks used inside doublediff()
p = params();
Hfloat = usurf / (1-p.rhoi/p.rhow);  % thickness if it were floating
float = (p.rhoi * Hfloat < p.rhow * (-topg));
icefree = (usurf < topg + 1.0);

fprintf('x range: %.3f km to %.3fkm\n',min(x)/1000,max(x)/1000)
fprintf('y range: %.3f km to %.3fkm\n',min(y)/1000,max(y)/1000)
fprintf('max thickness: %.3f m\n',max(max(thk)))
%fprintf('max thickness inconsistency: %.5f [frac]\n', max(max(abs(thk - (usurf - topg)))) / max(max(thk)))

Phi0arr = zeros(size(topg));
Phi = zeros(size(topg));
for i=1:length(x)
    for j=1:length(y)
        Phi0arr(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/p.spera; 
        if outline(i,j)>0.5, Phi(i,j) = Phi0arr(i,j); end
    end
end

if true
  fprintf('showing initial fields ...\n')

  sqw = 1000;  sqh = 310;  % sizes so pixels are square
  figure(1)
  set(gcf,'position',[100 250 sqw sqh])
  subplot(1,2,1)
  maxelev = max([max(max(topg)) max(max(usurf))]);
  minelev = min([min(min(topg)) min(min(usurf))]);
  imagesc(x/1000,flipud(-y)/1000,flipud(topg),[minelev maxelev]), colorbar
  title('bed elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(usurf),[minelev maxelev]), colorbar
  title('surface elevation  (m)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(2)
  set(gcf,'position',[100 200 sqw sqh])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(thk)), colorbar
  title('ice thickness  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  typemask = ones(size(thk));
  typemask(icefree) = 0.0;
  typemask(float) = 2.0;
  imagesc(x/1000,flipud(-y)/1000,flipud(typemask)), colorbar
  title('cell type:  0 = ice free,  1 = glacier,  2 = floating or ocean')
  xlabel('x (km)'), ylabel('y (km)')

  figure(3)
  set(gcf,'position',[100 150 sqw sqh])
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
[W, P] = doublediff(x,y,topg,usurf,magvb,outline,W0,P0,Phi,ts,te,5);

if true
  fprintf('showing final fields: W, P\n')

  figure(4)
  set(gcf,'position',[100 100 sqw sqh])
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

  figure(5)
  set(gcf,'position',[100 50 sqw sqh])
  subplot(1,2,1)
  imagesc(x/1000,flipud(-y)/1000,flipud(Po),[0 5.0e6]), colorbar
  title('overburden pressure P_o  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,flipud(-y)/1000,flipud(P),[0 5.0e6]), colorbar
  title('water pressure P  (Pa)')
  xlabel('x (km)'), ylabel('y (km)')

  figure(6)
  Wpos=W(W>0);  Ppos=P(W>0);  Popos=Po(W>0);
  plot(Wpos(:),Ppos(:)./Popos(:),'ko','markersize',8)
  WFC=0:.01:1.5*p.Wr;
  hold on,  plot(WFC,(WFC/p.Wr).^(7/2),'r--','linewidth',2.0), hold off
  %hold on,  plot(WFC,(WFC/(p.Wr/3)).^(7/2),'r--','linewidth',2.0), hold off
  axis([0 1.5*p.Wr 0 1.2])
  xlabel('water thickness W  (m)')
  ylabel('water pressure P as fraction of overburden P_o')
  %title('is P=P(W)?  (Flowers & Clarke is red dashed)')
  %print -dpdf isPofW.pdf
end

