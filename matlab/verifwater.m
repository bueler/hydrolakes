function [err, W, P] = verifwater(tyears,M,dofigs)
% VERIFWATER  Runs a verification test on a square region with  M  grid spaces
% in each direction for  tyears  years.  Initial condition is the result of
% RADIALSTEADY, which is a nearly-exact solution.  Thus measures drift away
% from a steady state.  Prints L_1 and L_inf errors for both W and P.
% Form:
%    [err, W, P] = verifwater(tyears,M,dofigs)
% where:
%    FIXME
% Calls:   RADIALSTEADY, DOUBLEDIFF

if nargin<1, tyears=1.0; end
if nargin<2, M=50; end
if nargin<3, dofigs=true; end

p = params();

% grid
Lx = 30.0e3;  % Lx > L in radialsteady()
dx = Lx / M;
dy = dx;
x = -Lx:dx:Lx;
y = x;

% dense radial grid and call to radialsteady()
fprintf('computing exact quantities as function of r with ODE solver ...\n')
fprintf('showing exact quantities as function of r ...\n')
%[r,Wrad,Prad,hrad,vbrad] = radialsteady(dofigs);  % use h0 and v0 defaults
[r,Wrad,Prad,hrad,vbrad] = radialsteady(false);  % use h0 and v0 defaults
R0 = 25.0e3;  % must be consistent with constant in RADIALSTEADY
L = 0.9 * R0;

% extend radial grid stuff so that piecewise-linear interpolation will work
r     = [50.0e3;  L+1;  r];
Wrad  = [0;       0;    Wrad];
Prad  = [0;       0;    Prad];
hrad  = [0;       0;    hrad];
vbrad = [0;       0;    vbrad];

% lookup table action:  get gridded data and exact quantities
[xx, yy] = ndgrid(x,y);
rr  = sqrt(xx.^2 + yy.^2);
h   = interp1(r,hrad,rr,'linear');
vb  = interp1(r,vbrad,rr,'linear');
WEX = interp1(r,Wrad,rr,'linear');
PEX = interp1(r,Prad,rr,'linear');
Po  = p.rhoi * p.g * h;
PEX = min(PEX,Po);

if dofigs
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
  imagesc(x/1000,y/1000,WEX,[0 p.Wr]), colorbar
  title('exact water thickness W  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,PEX / 1.0e5), colorbar
  title('exact water pressure P  (bar)')
  xlabel('x (km)'), ylabel('y (km)')
end

fprintf('calling doublediff() to do run for %.3f years on %d x %d grid ...\n',...
        tyears,M+1,M+1)
outline = ((abs(xx) < L) & (abs(yy) < L));
te = tyears * p.spera;
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a
Phi = Phi0 * ones(size(h));
PoL = p.rhoi * p.g * hrad(1);
Pfreebc = psteady(p,PoL,vbrad(1),Wrad(1));

% run silent:
[W, P] = doublediff(x,y,zeros(size(h)),h,vb,outline,WEX,PEX,Phi,0.0,te,5,true,Pfreebc);

if dofigs
  fprintf('showing numerical solution fields: W, P\n')

  figure(3)
  set(gcf,'position',[100 400 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,y/1000,W,[0 p.Wr]), colorbar
  title('numerical water thickness W  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,P / 1.0e5), colorbar
  title('numerical water pressure P  (bar)')
  xlabel('x (km)'), ylabel('y (km)')

  fprintf('showing error fields ...\n')

  figure(4)
  set(gcf,'position',[100 500 1200 400])
  subplot(1,2,1)
  imagesc(x/1000,y/1000,W-WEX), colorbar
  title('numerical water thickness error W - W_{exact}  (m)')
  xlabel('x (km)'), ylabel('y (km)')
  subplot(1,2,2)
  imagesc(x/1000,y/1000,(P-PEX)/1.0e5), colorbar
  title('numerical water pressure error P - P_{exact}  (bar)')
  xlabel('x (km)'), ylabel('y (km)')
end

Werrone = sum(sum(abs(W - WEX))) / (M+1)^2;
Perrone = sum(sum(abs(P - PEX))) / (M+1)^2;
Werrinf = max(max(abs(W - WEX)));
Perrinf = max(max(abs(P - PEX)));
WEXmax = max(max(WEX));  % not equal to  norm(WEX,'inf')
PEXmax = max(max(PEX));

fprintf('results:\n')
fprintf('          av  |W-Wexact| = %.6f m    av  |P-Pexact| = %.6f bar\n',...
        Werrone,Perrone/1e5)
fprintf('          [relative = %.6f]          [relative = %.6f]\n',...
        Werrone/WEXmax,Perrone/PEXmax)
fprintf('          max |W-Wexact| = %.6f m    max |P-Pexact| = %.6f bar\n',...
        Werrinf,Perrinf/1e5)
fprintf('          [relative = %.6f]          [relative = %.6f]\n',...
        Werrinf/WEXmax,Perrinf/PEXmax)

err = [Werrone Perrone Werrinf Perrinf];

