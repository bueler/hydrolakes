function [r,W] = radialsteady(M)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

p = params();

% constants specific to exact solution
h0   = 1000.0;         % m
R0   = 45.0e3;         % m
R1   = 20.0e3;         % m
Phi0 = 0.2 / p.spera;    % m/s
v0   = 100.0 / p.spera;  % m/s
vphi0 = Phi0 / (2 * p.c0);
CC   = (p.c1 * v0 / (p.c2 * p.A))^(1/3);

% radial grid
dr = R0 / M;
r = 0:dr:R0;

% ice thickness and sliding velocity
h = h0 * (1 - (r/R0).^2);
vb = v0 * (r - R1).^3 / (R0-R1)^3;
vb(r < R1) = 0.0;
figure(1), clf
ax = plotyy(r/1000.0,h,r/1000.0,vb * p.spera);
xlabel('r  (km)')
ylabel(ax(1),'h = ice thickness  (m)')
ylabel(ax(2),'|v_b| = ice sliding speed   (m/a)')

% solve the ODE; in octave this requires "odepkg"
wopt = odeset('RelTol', 1e-6,...
              'AbsTol', 1e-6);
%              'AbsTol', 1e-3,...
%              'InitialStep', -dr,...
%              'MaxStep', +dr);
%            "NormControl", "on", "OutputFcn", @odeplot);
[rr,WW] = ode45(@WODE,[R0 0.0],p.Wr-0.001,wopt);
fprintf('  [numerical ODE solution generated %d r values in 0 <= r <= R0]\n',length(rr))
figure(2), clf
plot(rr/1000.0,WW);
xlabel('r  (km)'), ylabel('W  (m)')

figure(3), clf
clear vb
vb = v0 * (rr - R1).^3 / (R0-R1)^3;
vb(rr < R1) = 0.0;
Po = p.rhoi * p.g * h0 * (1 - (rr/R0).^2);
P = Po - ( p.c1 * vb .* (p.Wr - WW) ./ (p.c2 * p.A * WW) ).^(1/3);
plot(rr/1000.0,Po/1e5,rr/1000.0,P/1e5);
xlabel('r  (km)'), ylabel('W  (m)')
legend('P_o = overburden pressure  (bar)','P = water pressure  (bar)')

% FIXME: extra stuff to look at return from ODE solver
%figure(99), clf
%rdens = zeros(1,length(r)-1);
%for j=1:length(r)-1,  rdens(j) = sum( ((r(j)<=rr) & (rr<r(j+1))) );  end
%plot(r(1:end-1)/1000.0,rdens,'.')
%xlabel('r  (km)'), ylabel('number of points near r')

  function dW = WODE(r,W)
  % this function defines the right side of the ODE for W
    if r < R1
      sb  = 0.0;
      dsb = 0.0;
    else
      dsb = CC / (R0 - R1);
      sb  = dsb * (r - R1);
    end
    dPodr = - (2 * p.rhoi * p.g * h0 / R0^2) * r;
    tmp1  = W.^(1/3) .* (p.Wr - W).^(2/3);
    numer = dsb * W * (p.Wr - W) - ( vphi0 * r + dPodr * W ) * tmp1;
    dW    = numer / ( (1/3) * p.Wr * sb + p.rhow * p.g * W * tmp1 );
  end

end

