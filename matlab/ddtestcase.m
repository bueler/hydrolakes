function ddtestcase(tyears,dx)
% DDTESTCASE  tests DOUBLEDIFF with artificial flat bed geometry

if nargin < 1
  tyears = 5.0;
end
if nargin < 2
  dx = 2.0e3;
end
if dx < 100.0, error('unreasonably fine grid'), end

% physical constants
spera = 31556926.0;
rhoi  = 910.0;         % kg m-3
rhow  = 1028.0;        % kg m-3
g     = 9.81;          % m s-2

% square of side length 100 km
L = 100.0e3;
x = 0:dx:L;
y = 0:dx:L;
[xx,yy] = ndgrid(x,y);

fprintf('running DDTESTCASE on %d x %d grid with spacing dx = %.3f km ...\n',...
        length(x),length(y),dx/1000)

% build geometry from 1/4 of Halfar cap at t0; see Karthaus notes
H0 = 1000.0;   % m
R0 = 90.0e3;  % = 90 km
A = 1.0e-16 / spera;
Gamma  = 2 * A * (rhoi * g)^3 / 5;
t0 = (18*Gamma)^(-1) * (7/4)^3 * (R0^4/H0^7);
fprintf('Halfar geometry with t0 = %.4e a\n',t0/spera)
rr = sqrt(xx.*xx + yy.*yy);
rr = rr / R0;
h = H0 * max(0, 1 - rr.^(4/3)).^(3/7);
b = zeros(size(h));

% water source is 1 m/a in patch around (xw,yw)
xw = 0.45 * R0;  yw = xw;
Rw = 0.3 * R0;
wpatch = ( (xx - xw).^2 + (yy - yw).^2 < (0.3 * R0)^2 );
Phi = (1.0/spera) * wpatch;

% set up other fields
Po = rhoi * g * h;
W0 = zeros(size(h));
magvb = zeros(size(h));
outline = (h > 0);

% initial fields: surface elevation and Phi
figure(91)
imagesc(x/1000,y/1000,Po / 1e5), colorbar
title('DDTESTCASE: overburden pressure = initial pressure  (bar)')
xlabel('x (km)'), ylabel('y (km)')

figure(92)
imagesc(x/1000,y/1000,Phi*spera), colorbar
title('DDTESTCASE: water input (m/a)')
xlabel('x (km)'), ylabel('y (km)')

% run doublediff
fprintf('running DOUBLEDIFF ...\n')
%[W, P] = doublediff(x,y,b,h,magvb,outline,W0,Po,Phi,0.0,tyears*spera,10);
[W, P] = doublediff(x,y,b,h,magvb,outline,W0,0.8*Po,Phi,0.0,tyears*spera,10);

% final fields: water thickness and pressure
figure(93)
imagesc(x/1000,y/1000,W), colorbar
title('DDTESTCASE: final water thickness (m)')
xlabel('x (km)'), ylabel('y (km)')

figure(94)
imagesc(x/1000,y/1000,P / 1e5), colorbar
title('DDTESTCASE: final pressure  (bar)')
xlabel('x (km)'), ylabel('y (km)')

figure(95)
imagesc(x/1000,y/1000,(Po - P) / 1e5), colorbar
title('DDTESTCASE: final effective pressure  (bar)')
xlabel('x (km)'), ylabel('y (km)')

