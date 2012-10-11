function [r,W] = radialsteady(M)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

spera = 31556926.0;
rhoi  = 910.0;         % kg m-3
rhow  = 1028.0;        % kg m-3
g     = 9.81;          % m s-2

% major model parameters:
A  = 3.1689e-24;       % ice softness (Pa-3 s-1)
K  = 1.0e-2;           % m s-1   FIXME: want Kmax or Kmin according to W > Wr
Wr = 1.0;              % m
c1 = 0.500;            % m-1
c2 = 0.040;            % [pure]
c0 = K / (rhow * g);

% constants specific to exact solution
h0   = 1000.0;         % m
R0   = 45.0e3;         % m
R1   = 20.0e3;         % m
Phi0 = 0.1 / spera;    % m/s
v0   = 100.0 / spera;  % m/s
varphi0 = Phi0 / (2 * c0);
CC   = (c1 * v0 / (c2 * A))^(1/3);

% radial grid
dr = R0 / M;
r = 0:dr:R0;

% ice thickness and sliding velocity
h = h0 * (1 - (r/R0).^2);
vb = v0 * (r - R1).^3 / (R0-R1)^3;
vb(r < R1) = 0.0;
figure(1), clf
ax = plotyy(r/1000.0,h,r/1000.0,vb * spera);
xlabel('r  (km)')
ylabel(ax(1),'h = ice thickness  (m)')
ylabel(ax(2),'|v_b| = ice sliding speed   (m/a)')

% solve the ODE; in octave this requires "odepkg"
wopt = odeset('RelTol', 1e-3,...
              'AbsTol', 1e-3,...
              'InitialStep', -dr,...
              'MaxStep', +dr);
       %            "NormControl", "on", "OutputFcn", @odeplot);
[rr,WW] = ode45(@WODE,[R0 0],Wr,wopt);
size(rr), size(WW)
figure(1), clf
plot(rr/1000.0,WW);
xlabel('r  (km)')
rr
WW

  function dW = WODE(rr,W)
  % this function defines the right side of the ODE for W
    if rr < R1
      sb  = 0.0;
      dsb = 0.0;
    else
      dsb = CC / (R0 - R1);
      sb  = dsb * (rr - R0);
    end
    dPodr = - (2 * rhoi * g * h0 / R0^2) * rr;
    tmp1  = W.^(1/3) .* (Wr - W).^(2/3);
    tmp2  = dsb * W * (Wr - W) - ( varphi0 + dPodr * W ) * tmp1;
    dW    = tmp2 / ( (1/3) * Wr * sb + rhow * g * W * tmp1 );
  end

