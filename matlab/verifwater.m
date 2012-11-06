function [W, P] = verifwater(tyears,M)
% VERIFWATER  **FIXME this comment**
% Form:  [W, P] = verifwater(tyears,M)
% Calls:   RADIALSTEADY, DOUBLEDIFF

if nargin<1, tyears=5.0; end
if nargin<2, M=50; end

p = params();

% grid
L = 30.0e3;  % L > R0 in radialsteady()
dx = L / M;
dy = dx;
x = -L:dx:L;
y = x;

% dense radial grid and call to radialsteady()
[r,Wrad,hrad,vbrad] = radialsteady(M*10);  % use h0 and v0 defaults
R0 = 25.0e3;  % must be consistent
r    = [r     (R0+10) 50.0e3];
Wrad = [Wrad  0       0];
hrad = [hrad  0       0];
vb   = [vbrad 0       0];

% lookup
[xx, yy] = ndgrid(x,y);
rr  = sqrt(xx.^2 + yy.^2);
thk = interp1(r,hrad,rr,'linear');
vb  = interp1(r,vbrad,rr,'linear');

if true
  fprintf('showing initial fields ...\n')

  figure(1)
  set(gcf,'position',[100 200 1200 400])
  subplot(1,2,1)
  imagesc(xx/1000,yy/1000,thk), colorbar
  title('ice thickness  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(xx/1000,yy/1000,vb), colorbar
  title('ice sliding velocity')
  xlabel('x (km)'), ylabel('y (km)')
end

% FIXME: proceed to call doublediff()


