function [r,W,P,h,vb] = radialsteady(dofigs)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

if nargin<1, dofigs=true; end

p = params();

% constants specific to exact solution
h0   = 500.0;           % m     center thickness
v0   = 100.0 / p.spera; % m/s   sliding velocity at margin
R0   = 25.0e3;           % m     ice sheet margin
R1   = 5.0e3;            % m     onset of sliding
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a

vphi0 = Phi0 / (2 * p.c0);

if dofigs
  % ice thickness and sliding velocity for plot
  dr = R0 / 500;
  rr = 0:dr:R0;
  hrr = h0 * (1 - (rr/R0).^2);
  vbrr = v0 * (rr - R1).^5 / (R0-R1)^5;
  vbrr(rr < R1) = 0.0;

  set(0,'defaultaxesfontsize',14)
  set(0,'defaultlinelinewidth',3.0)
  figure(97), clf
  ax = plotyy(rr/1000.0,hrr,rr/1000.0,vbrr * p.spera);
  xlabel('r  (km)')
  ylabel(ax(1),'h = ice thickness  (m)')
  ylabel(ax(2),'|v_b| = ice sliding speed   (m/a)')

  clear dr rr hrr vbrr
end

% in octave this requires "odepkg", but then fails because of negative step direction

% solve the ODE
wopt = odeset('RelTol', 1e-12,'AbsTol', 1e-8);
% others: 'InitialStep', 'MaxStep', 'NormControl', 'OutputFcn'
[r,W] = ode45(@WODE,[R0 0.0],p.Wr-0.001,wopt);
fprintf('  [numerical ODE solution generated %d r values in 0 <= r <= R0]\n',length(r))

if dofigs
  % show ODE soln W(r)
  figure(98), clf
  plot(r/1000.0,W);
  xlabel('r  (km)'), ylabel('W  (m)')
  axis([0 R0/1000.0 0 1.2*p.Wr]), grid on
  title(sprintf('exact solution W(r)   (note W_r = %.2f m)\n',p.Wr))
end

% compute pressure solution
vb = v0 * (r - R1).^5 / (R0-R1)^5;
vb(r < R1) = 0.0;
vb(r > R0) = 0.0;
h = h0 * (1 - (r/R0).^2);
h(r > R0) = 0.0;
Po = p.rhoi * p.g * h;
P = Psteady(Po,vb,W);

if dofigs
  % show pressure solution
  figure(99), clf
  plot(r/1000.0,Po/1e5,r/1000.0,P/1e5);
  xlabel('r  (km)'), ylabel('pressure  (bar)')
  legend('P_o(r) = overburden pressure','P(r) = exact water pressure')
end

  function dW = WODE(r,W)
  % this function defines the right side of the ODE for W:  W'(r) = WODE(r,W)
  if r < R1
    sb  = 0.0;
    dsb = 0.0;
  else
    CC  = p.c1 / (p.c2 * p.A);
    % vb = v0 * (r - R1).^5 / (R0-R1)^5   and   sb = (CC * vb)^(1/3)
    CZ  = (CC * v0)^(1/3);
    zz  = ((r - R1) / (R0 - R1)).^(1/3);
    sb  = CZ * zz.^5;
    CD  = (5 * CZ) / (3 * (R0 - R1));
    dsb = CD * zz.^2;
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

