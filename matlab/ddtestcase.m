% DDTESTCASE  tests DOUBLEDIFF with artificial flat bed geometry

% physical constants
spera = 31556926.0;
rhoi  = 910.0;         % kg m-3
rhow  = 1028.0;        % kg m-3
g     = 9.81;          % m s-2

% square of side length 100 km
L = 100.0e3;
dx = 2.0;
x = 0:dx:L;
y = 0:dy:L;
[xx,yy] = meshgrid(x,y);

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

% set up fields
Po = rhoi * g * h;
W0 = zeros(size(h));
magvb = zeros(size(h));
outline = (h > 0);
Phi = FIXME

% run doublediff
[W, P] = doublediff(x,y,b,h,magvb,outline,W0,Po,Phi,0.0,1.0*spera,10)

