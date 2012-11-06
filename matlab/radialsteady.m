function [r,W,h,vb] = radialsteady(M,h0,v0)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

p = params();

% defaults
if nargin < 1
  M = 100;
end
if nargin < 2
  h0   = 500.0;           % m     center thickness
end
if nargin < 3
  v0   = 100.0 / p.spera; % m/s   sliding velocity at margin
end

% constants specific to exact solution
R0   = 25.0e3;           % m     ice sheet margin
R1   = 5.0e3;            % m     onset of sliding
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a

vphi0 = Phi0 / (2 * p.c0);

% radial grid
dr = R0 / M;
r = 0:dr:R0;

% ice thickness and sliding velocity
h = h0 * (1 - (r/R0).^2);
vb = v0 * (r - R1).^4 / (R0-R1)^4;
vb(r < R1) = 0.0;

% show
figure(97), clf
ax = plotyy(r/1000.0,h,r/1000.0,vb * p.spera);
xlabel('r  (km)')
ylabel(ax(1),'h = ice thickness  (m)')
ylabel(ax(2),'|v_b| = ice sliding speed   (m/a)')

% in octave this requires "odepkg", but then fails because of negative step direction

% solve the ODE
wopt = odeset('RelTol', 1e-6,'AbsTol', 1e-6);
% others: 'InitialStep', 'MaxStep', 'NormControl', 'OutputFcn'
[rr,WW] = ode45(@WODE,[R0 0.0],p.Wr-0.001,wopt);
fprintf('  [numerical ODE solution generated %d r values in 0 <= r <= R0]\n',length(rr))

% show ODE soln W(r)
figure(98), clf
plot(rr/1000.0,WW);
xlabel('r  (km)'), ylabel('W  (m)')

% show pressure solution
figure(99), clf
vbrr = v0 * (rr - R1).^4 / (R0-R1)^4;  % vb calculation using ODE rr values
vbrr(rr < R1) = 0.0;
Porr = p.rhoi * p.g * h0 * (1 - (rr/R0).^2);
plot(rr/1000.0,Porr/1e5,rr/1000.0,Psteady(Porr,vbrr,WW)/1e5);
xlabel('r  (km)'), ylabel('pressure  (bar)')
legend('P_o = overburden pressure','P = water pressure')

  function dW = WODE(r,W)
  % this function defines the right side of the ODE for W:  W'(r) = WODE(r,W)
  if r < R1
    sb  = 0.0;
    dsb = 0.0;
  else
    CC  = p.c1 / (p.c2 * p.A);
    % vb = v0 * (r - R1).^4 / (R0-R1)^4   and   sb = (CC * vb)^(1/3)
    CZ  = (CC * v0)^(1/3);
    zz  = ((r - R1) / (R0 - R1)).^(1/3);
    sb  = CZ * zz.^4;
    CD  = (4 * CZ) / (3 * (R0 - R1));
    dsb = CD * zz;
  end
  dPodr = - (2 * p.rhoi * p.g * h0 / R0^2) * r;
  tmp1  = W.^(1/3) .* (p.Wr - W).^(2/3);
  numer = dsb .* W .* (p.Wr - W) - ( vphi0 * r + dPodr .* W ) .* tmp1;
  denom = (1/3) * p.Wr * sb + p.rhow * p.g * W .* tmp1;
  dW    = numer ./ denom;
  end

  function P = Psteady(Po,vb,W)
  % this function computes P(W) in steady state
  if any(Po < 0), error('Psteady() requires nonnegative overburden pressure Po'), end
  if any(vb < 0), error('Psteady() requires nonnegative sliding speed vb'), end
  if any(W <= 0), error('Psteady() requires positive water thickness W'), end
  ff = (p.Wr - W) ./ W;
  ff(ff < 0) = 0.0;
  CC = p.c1 / (p.c2 * p.A);
  P  = Po - (CC * vb .* ff).^(1/3);
  P(P < 0) = 0.0;
  end

end

