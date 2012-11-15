function [r,W,P,h,vb] = radialsteady(dofigs)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

if nargin<1, dofigs=true; end

p = params();

% constants specific to exact solution
h0   = 500.0;            % m     center thickness
v0   = 100.0 / p.spera;  % m/s   sliding velocity at margin
R0   = 25.0e3;           % m     ice sheet margin
R1   = 5.0e3;            % m     onset of sliding
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a

vphi0 = Phi0 / (2 * p.c0);
L = 0.9 * R0;

% WcL is key constant in construction
hL  = h0 * (1 - (L/R0).^2);
vbL = v0 * (L - R1).^5 / (R0-R1)^5;
PoL = p.rhoi * p.g * hL;
sbL = ( (p.c1 * vbL / (p.c2 * p.A)) )^(1/3);
WcL = (sbL^3 * p.Wr - PoL^3 * p.Y0) / (sbL^3 + PoL^3)

if dofigs
  set(0,'defaultaxesfontsize',14)
  set(0,'defaultlinelinewidth',3.0)

  % ice thickness and sliding velocity for plot
  dr = R0 / 500;  r = 0:dr:R0;
  h = h0 * (1 - (r/R0).^2);
  vb = v0 * (r - R1).^5 / (R0-R1)^5;
  vb(r < R1) = 0.0;
  figure(96), clf
  ax = plotyy(r/1000.0,h,r/1000.0,vb * p.spera);
  xlabel('r  (km)')
  ylabel(ax(1),'h = ice thickness  (m)')
  ylabel(ax(2),'|v_b| = ice sliding speed   (m/a)')

  % grid for showing O, N, U regions behind exact W(r) solution
  r = 0:dr:L;
  dW = p.Wr/500;
  W = dW:dW:1.2*p.Wr;
  [rr WW] = meshgrid(r,W);
  Po = p.rhoi * p.g * h0 * (1 - (rr/R0).^2);
  vb = v0 * (rr - R1).^5 / (R0-R1)^5;
  vb(rr < R1) = 0.0;
  ratio = Psteady(p,Po,vb,WW) ./ Po;
  figure(97), clf
  [cc ll] = contour(rr/1000.0,WW,ratio,[0.000001 0.999999],'k','linewidth',1.5);
  %clabel(cc,ll,'fontsize',14.0);
  text(2,0.5,'O','fontsize',16.0)
  text(12,1.1,'O','fontsize',16.0)
  text(12,0.5,'N','fontsize',16.0)
  text(22,0.3,'U','fontsize',16.0)
  xlabel('r  (km)'), ylabel('W  (m)')
  axis([0 R0/1000.0 0 1.2*p.Wr]), grid on

  clear dr r rr h Po vb W WW
end

return
% in octave this requires "odepkg", but then fails because of negative step direction

% solve the ODE
wopt = odeset('RelTol', 1e-12,'AbsTol', 1e-8);
% others: 'InitialStep', 'MaxStep', 'NormControl', 'OutputFcn'
[r,W] = ode45(@WODE,[R0 0.0],p.Wr-0.001,wopt);
fprintf('  [numerical ODE solution generated %d r values in 0 <= r <= R0]\n',length(r))

if dofigs
  % show ODE soln W(r) onto existing fig
  figure(97), hold on
  plot(r/1000.0,W,'k--');
  hold off
end

% compute pressure solution
vb = v0 * (r - R1).^5 / (R0-R1)^5;
vb(r < R1) = 0.0;
vb(r > R0) = 0.0;
h = h0 * (1 - (r/R0).^2);
h(r > R0) = 0.0;
Po = p.rhoi * p.g * h;
P = Psteady(p,Po,vb,W);

if dofigs
  % show pressure solution
  figure(98), clf
  plot(r/1000.0,Po/1e5,'k',r/1000.0,P/1e5,'k--');
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

  function P = Psteady(p,Po,vb,W)
  % this function computes P(W) in steady state
  % it is vectorized in all arguments
  if any(any(Po < 0))
    error('Psteady() requires nonnegative overburden pressure Po'), end
  if any(any(vb < 0))
    error('Psteady() requires nonnegative sliding speed vb'), end
  if any(any(W <= 0))
    error('Psteady() requires positive water thickness W'), end
  sbcube = p.c1 * vb / (p.c2 * p.A);
  frac = max(0.0, p.Wr - W) ./ (W + p.Y0);
  P = Po - (sbcube .* frac).^(1/3);
  P(P < 0) = 0.0;
  end

end

