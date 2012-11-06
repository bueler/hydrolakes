function [W, P] = verifwater(tyears,M)
% VERIFWATER  Runs a verification test on a square region with  M  grid spaces
% in each direction for  tyears  years.  Initial condition is the result of
% RADIALSTEADY, which is a nearly-exact solution.  Thus measures drift away
% from a steady state.  Prints L_1 and L_inf errors for both W and P.
% Form:  [W, P] = verifwater(tyears,M)
% Calls:   RADIALSTEADY, DOUBLEDIFF

if nargin<1, tyears=1.0; end
if nargin<2, M=50; end

p = params();

% grid
L = 30.0e3;  % L > R0 in radialsteady()
dx = L / M;
dy = dx;
x = -L:dx:L;
y = x;

% dense radial grid and call to radialsteady()
fprintf('computing exact quantities as function of r with ODE solver ...\n')
fprintf('showing exact quantities as function of r ...\n')
[r,Wrad,Prad,hrad,vbrad] = radialsteady();  % use h0 and v0 defaults
R0 = 25.0e3;  % must be consistent with constant in RADIALSTEADY

% extend radial grid stuff so that piecewise-linear interpolation will work
r     = [r;     (R0+1); 50.0e3];
Wrad  = [Wrad;  0;       0];
Prad  = [Prad;  0;       0];
hrad  = [hrad;  0;       0];
vbrad = [vbrad; 0;       0];

% lookup table action:  get gridded data and exact quantities
[xx, yy] = ndgrid(x,y);
rr  = sqrt(xx.^2 + yy.^2);
h   = interp1(r,hrad,rr,'linear');
vb  = interp1(r,vbrad,rr,'linear');
WEX = interp1(r,Wrad,rr,'linear');
PEX = interp1(r,Prad,rr,'linear');
Po  = p.rhoi * p.g * h;
PEX = min(PEX,Po);

if true
  fprintf('showing data fields: h, |v_b|\n')

  figure(1)
  set(gcf,'position',[100 200 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,y/1000,h), colorbar
  title('ice thickness h (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,vb * p.spera), colorbar
  title('ice sliding speed vb (m/a)')
  xlabel('x (km)'), ylabel('y (km)')

  fprintf('showing exact solution fields: W, P\n')

  figure(2)
  set(gcf,'position',[100 300 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,y/1000,WEX), colorbar
  title('exact water thickness W  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,PEX / 1.0e5), colorbar
  title('exact water pressure P  (bar)')
  xlabel('x (km)'), ylabel('y (km)')
end

fprintf('calling doublediff() to do run for %.3f years ...\n\n',tyears)
outline = ((abs(xx) < L) & (abs(yy) < L));
te = tyears * p.spera;
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a
Phi = Phi0 * ones(size(h));
%[W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin)
[W, P] = doublediff(x,y,zeros(size(h)),h,vb,outline,WEX,PEX,Phi,0.0,te,5);

if true
  fprintf('showing error fields ...\n')

  figure(3)
  set(gcf,'position',[100 300 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,y/1000,W-WEX), colorbar
  title('numerical water thickness error W - W_{exact}  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,(P-PEX)/1.0e5), colorbar
  title('numerical water pressure error P - P_{exact}  (bar)')
  xlabel('x (km)'), ylabel('y (km)')
end

