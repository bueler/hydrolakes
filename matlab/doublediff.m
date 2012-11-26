function [W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin,silent)
% DOUBLEDIFF  A nearly-full-cavity, Darcy flux subglacial hydrology model.
% This model combines an advection-diffusion equation for water thickness W
% with a diffusion equation for the pressure P.
%
% MODEL EQUATIONS:  The equations are described in dampnotes.pdf.
% Here is a summary only.  The state space functions are W and P.
% Derived fields are:
%   psi = P + rhow g b                hydraulic potential
%   V   = - c0 grad P - K grad b      velocity of water
%   P_o = rhoi g H                    overburden pressure
% Note that the velocity V has components  (alpha,beta)  in the code.
% Evolution equations are
%   dW/dt + div( V W ) = div ( K W grad W ) + Phi
%   (E0/P_o) dP/dt = div ( c0 W grad psi ) + Clos(P,W) - Open(W) + Phi
% where
%   Open(W) = c1 magvb (Wr - W)_+
%   Clos(P,W) = c2 A (Po - P)^3 (W + Y0)
% Various constants come from PARAMS.
%
% USAGE:  The calling form is
%
%      [W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin,silent)
%
% The input data are given on (Mx+1) x (My+1) grids:
%   x = coordinate vector of length Mx+1
%   y = coordinate vector of length My+1
%   b = bed elevation; array of (Mx+1) x (My+1) values
%   h = surface elevation; <same size>
%   magvb = |v_b| = magnitude of sliding velocity; <same size>
%   outline = 1 if inside modeled region, 0 otherwise; <same size>
%   W0 = initial values for water thickness W in m; <same size>
%   P0 = initial values for pressure P in Pa; <same size>
%   Phi = melt rate as m s-1; <same size>
%
% Other parameters last:
%   ts = run start time
%   te = run end time
%   Nmin = use at least this many time steps.
%   silent = if true then suppress all stdout
%
% The maximum time step is  (te-ts)/Nmin,  but the code does adaptive
% time-stepping.  At each time step the code reports time step, time step
% restrictions (CFL for advection and criterion for W diffusion and
% P diffusion), and total water amount.
%
% Note that where outline==0 we set the water thickness W to zero.
%
% The output variables are the final values of the state variables W,P.
%
% See also PARAMS, PLOTPSTEADY, RADIALSTEADY, VERIFWATER, NBREENWATER.

p = params();
if nargin<13, silent=false; end

% get grid parameters from W0; other fields must match but code does not check
[Mx, My] = size(W0);   Mx = Mx-1;   My = My-1;
dx = x(2) - x(1);
dy = y(2) - y(1);

% generate staggered grid
xs = 0.5 * (x(1:end-1) + x(2:end));
ys = 0.5 * (y(1:end-1) + y(2:end));

% bed slope components onto staggered grid
dbdx = (b(2:end,:) - b(1:end-1,:)) / dx;
dbdy = (b(:,2:end) - b(:,1:end-1)) / dy;

% initialize timestepping
dtmax = (te - ts) / Nmin;
t = ts;

% initialize state variable W and mass accounting
if any(any(W0<0))
  error('points exist where initial water thickness is negative'), end
W = W0;
volW = 0.0;
dA = dx*dy;

% overburden pressure will not vary during run
% note b is bedrock elevation, not ice lower surface (where floating)
if any(any(h < b))
  error('points exist where ice surface is below bed'), end
if any(any(h < 0))
  error('points exist where ice surface is below sea level'), end
Po = p.rhoi * p.g * (h - b);
rhosw = 1028.0;
Hfloat = h / (1 - p.rhoi / rhosw);  % h=surfelev; Hfloat = thickness if it were floating
float = (p.rhoi * Hfloat < rhosw * (-b));
Po(float) = p.rhoi * p.g * Hfloat(float);  % if floating, reset to ice base pressure

% is initial state valid?
if any(any(P0<0))
  error('points exist where initial pressure is negative'), end
if any(any(P0>Po))
  error('points exist where initial pressure is above overburden'), end

% impose P boundary conditions based on geometry only and initialize P
P0(float) = Po(float);
icefree = (h < b + 1.0) & (h >= 0.0);
P0(icefree) = 0.0;
P = P0;

% locations outside subglacial aquifer where W and P are known
known = (icefree | float);

if ~silent
  fprintf('running ...\n\n')
  fprintf('   [max V (m/a)   dtCFL (a)   dtDIFF_W (a)   dtDIFF_P (a)  -->  dt (a)]\n')
  fprintf('t (a):   max W (m)  av W (m)  vol(km^3)  rel-bal\n\n')
end

while t<te
  % grad pressure
  dPdx = (P(2:end,:) - P(1:end-1,:)) / dx;
  dPdy = (P(:,2:end) - P(:,1:end-1)) / dy;

  % velocity  V = - c0 grad P - K grad b = (alphV,betaV)
  alphV = - p.c0 * dPdx - p.K * dbdx;
  betaV = - p.c0 * dPdy - p.K * dbdy;

  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));

  % determine time step adaptively
  % FIXME: dtCFL can be infinity if velocity is zero because P and b are constant
  dtCFL = 0.5 / (max(max(abs(alphV)))/dx + max(max(abs(betaV)))/dy);
  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFFW = 0.25 / (p.K * maxW * (1/dx^2 + 1/dy^2));
  maxH = max(max(h-b)) + 2 * p.E0;  % regularized: forces dtDIFFP < dtDIFFW
  dtDIFFP = (p.rhow * p.E0 / (p.rhoi * maxH)) * dtDIFFW;
  dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);

  % report on time step
  maxV = sqrt(max(max(alphV))^2 + max(max(betaV))^2);
  if ~silent
    fprintf('   [%.5e  %.6f  %.6f  %.6f   -->  dt = %.6f (a)]\n',...
            maxV*p.spera, dtCFL/p.spera, dtDIFFW/p.spera, dtDIFFP/p.spera, dt/p.spera)
  end

  % hydraulic potential
  psi = P + p.rhow * p.g * (b + W);
  psi(float) = Po(float);

  % opening and closure terms in pressure equation
  Open = p.c1 * magvb .* (p.Wr - W);
  Open(Open < 0.0) = 0.0;
  Clos = p.c2 * p.A * (Po - P).^3 .* (W + p.Y0);

  % P time step
  pux = p.c0 / dx^2;
  puy = p.c0 / dy^2;
  for i=2:length(x)-1
    for j=2:length(y)-1
      psiij = psi(i,j);
      tmp = 0;
      if ~known(i+1,j) & ~known(i-1,j)
        tmp = tmp + pux * (Wea(i,j) * (psi(i+1,j)-psiij) - Wea(i-1,j) * (psiij-psi(i-1,j)));
      end
      if ~known(i,j+1) & ~known(i,j-1)
        tmp = tmp + puy * (Wno(i,j) * (psi(i,j+1)-psiij) - Wno(i,j-1) * (psiij-psi(i,j-1)));
      end
      Ptmp = P(i,j) + (dt * Po(i,j) / p.E0) * ( tmp + Clos(i,j) - Open(i,j) + Phi(i,j) );
      % projection:
      Ptmp = max(0.0, Ptmp);
      Ptmp = min(Ptmp, Po(i,j));
      P(i,j) = Ptmp;
    end
  end
  P(float) = Po(float);
  P(icefree) = 0.0;

% DECOUPLE POINT
%P = P0;

  % idea: at this point we could recompute V=(alph,beta), but then dtCFL would be invalid

  % coefficients for W time step
  nux = dt / dx;
  nuy = dt / dy;
  mux = p.K * dt / dx^2;
  muy = p.K * dt / dy^2;

  % build Qe using either first-order or limiter
  Qe = zeros(size(alphV));
  for i=1:length(x)-1
    for j=1:length(y)
      ve = alphV(i,j);
      if ~limiter
        psi = 0.0;
      else
        if ve >= 0
          if (i == 1) | (W(i+1,j) == W(i,j))
            psi = 0;
          else
            theta = ( W(i,j) - W(i-1,j) ) / ( W(i+1,j) - W(i,j) );
            psi = psikoren(theta);
          end
        else
          if (i == length(x)-1) | (W(i+1,j) == W(i,j))
            psi = 0;
          else
            FIXME: check this!!
            thetainv = ( W(i+2,j) - W(i+1,j) ) / ( W(i+1,j) - W(i,j) );
            psi = psikoren(thetainv);
          end
        end
      end
      % now compute Q
      if ve >= 0
        Qe(i,j) = ve * ( W(i,j)   + psi * ( W(i+1,j) -   W(i,j) ) );
      else
        Qe(i,j) = ve * ( W(i+1,j) + psi * ( W(i,j)   - W(i+1,j) ) );
      end
    end
  end

  % build Qe using either first-order or limiter
  FIXME: to do
  Qn = zeros(size(betaV));
  for i=1:length(x)
    for j=1:length(y)-1
      vn = betaV(i,j);
      if limiter
      else
        if vn >= 0
          Qn(i,j) = vn * W(i,j);
        else
          Qn(i,j) = vn * W(i,j+1);
        end
      end
    end
  end

  % W time step
  Wnew = W;  % copies unaltered initial b.c.s
  inputvol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      Wij = W(i,j);
      tmp = 0;
      if ~known(i,j) & ~known(i+1,j)
        tmp = tmp + mux * Wea(i,j)   * (W(i+1,j)-Wij);
      end
      if ~known(i,j) & ~known(i-1,j)
        tmp = tmp - mux * Wea(i-1,j) * (Wij-W(i-1,j));
      end
      if ~known(i,j) & ~known(i,j+1)
        tmp = tmp + muy * Wno(i,j)   * (W(i,j+1)-Wij);
      end
      if ~known(i,j) & ~known(i,j-1)
        tmp = tmp - muy * Wno(i,j-1) * (Wij-W(i,j-1));
      end
      inputdepth = dt * Phi(i,j);
      Wnew(i,j) = Wij - nux * (Qe(i,j) - Qe(i-1,j)) - nuy * (Qn(i,j) - Qn(i,j-1)) ...
                    + tmp + inputdepth;
      inputvol = inputvol + inputdepth * dA;
    end
  end

  %FIXME: W < 0 is not an error here if Phi < 0, but needs fix: Wnew(i,j)=0
  if any(any(Wnew < 0)), error('negative W'), end

  % do water removal at obvious cells (and volume accounting)
  losevol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      if (outline(i,j) < 0.5) | known(i,j)
        losevol = losevol + Wnew(i,j)*dA;
        Wnew(i,j) = 0.0;
      end
    end
  end

  % report on new state W
  sumnew = sum(sum(Wnew));
  volnew = sumnew * dA;
  vbalance = (volW + inputvol - losevol) - volnew;
  if volnew > 0, relbal = abs(vbalance/volnew); else relbal = nan; end
  if ~silent
    fprintf('t = %.5f (a):  %7.3f  %7.3f        %.3f        %.0e\n',...
            (t+dt)/p.spera, max(max(Wnew)), sumnew/((Mx+1)*(My+1)),...
            volnew/(1e9),relbal)
  end

  % actually update to new state W
  volW = volnew;
  W = Wnew;

% DECOUPLE POINT
%W = W0;

  t = t + dt;
end

  function z = psikoren(theta)
    z = 0.0;  %FIXME
  end

