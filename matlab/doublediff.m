function [W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin)
% DOUBLEDIFF  A nearly-full-cavity, Darcy flux subglacial hydrology model.
% This model combines an advection-diffusion equation for water thickness W
% with a penalized form of the "Y=W" full-cavity condition.  The latter
% evolves the hydraulic potential psi by a diffusion equation.
%
% ** MODEL EQUATIONS **:  The equations are described in dampnotes.pdf.
% Here is a summary only.  The major ("state space") functions are W and P.
% Derived fields are:
%   V = - (K / (rhow g)) grad P - K grad b   % velocity of water
%   P_o = rhoi g H                           % overburden pressure
% Note that the velocity V has components  (alpha,beta)  in the code.
% The evolution equations are
%
%   dW/dt + div( V W ) = div ( K W grad W ) + Phi
%   (E0/P_o) dP/dt = FIXME
%
% ** USAGE **  The calling form is
%
%      [W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin)
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
% Regarding the last three scalar arguments, the code runs from times
% ts  to  te  using at least  Nmin  time steps.  Thus the maximum time step
% is  (te-ts)/Nmin.  But the code does adaptive time-stepping.  At each time
% step the code reports time step, time step restrictions (CFL for advection
% and criterion for diffusion), and total water amount.
%
% Note that where outline==0 we set the water thickness W to zero.
%
% The output variables are the final values of the state variables W,P.

p = params();

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
Hfloat = h / (1-p.rhoi/p.rhow);  % thickness if it were floating
float = (p.rhoi * Hfloat < p.rhow * (-b));
Po(float) = p.rhoi * p.g * Hfloat(float);  % if floating, set to ice base pressure

% initialize state variable P
if any(any(P0<0))
  error('points exist where initial pressure is negative'), end
if any(any(P0>Po))
  error('points exist where initial pressure is above overburden'), end
P = P0;

% impose P boundary conditions based on geometry only
P(float) = Po(float);
icefree = (h < b + 1.0);
P(icefree) = 0.0;

fprintf('running ...\n\n')
fprintf('   [max V (m/a)   dtCFL (a)   dtDIFF_W (a)   dtDIFF_P (a)  -->  dt (a)]\n')
fprintf('t (a):   max W (m)  av W (m)  vol(km^3)  rel-bal\n\n')

while t<te
  % hydraulic potential and terms in pressure equation
  psi = P + p.rhow * p.g * (b + W);
  Ocav = p.c1 * magvb .* (p.Wr - W);
  Ocav(Ocav < 0.0) = 0.0;
  Ccrp = p.c2 * p.A * (Po - P).^3 .* W;

  % grad pressure and grad water
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
  fprintf('   [%.5e  %.6f  %.6f  %.6f   -->  dt = %.6f (a)]\n',...
          maxV*p.spera, dtCFL/p.spera, dtDIFFW/p.spera, dtDIFFP/p.spera, dt/p.spera)

  % coefficients for time step
  nux = dt / dx;
  nuy = dt / dy;
  mux = p.K * dt / dx^2;
  muy = p.K * dt / dy^2;
  pux = p.c0 * dt / dx^2;
  puy = p.c0 * dt / dy^2;

  % P time step
  Pnew = P;   % will keep P at edge of computational domain unaltered
  for i=2:length(x)-1
    for j=2:length(y)-1
      psiij = psi(i,j);
      tmp = pux * (Wea(i,j) * (psi(i+1,j)-psiij) - Wea(i-1,j) * (psiij-psi(i-1,j))) + ...
            puy * (Wno(i,j) * (psi(i,j+1)-psiij) - Wno(i,j-1) * (psiij-psi(i,j-1)));
      Ptmp = P(i,j) + (Po(i,j) / p.E0) * ( tmp + dt * Ccrp(i,j) - ...
                                         dt * Ocav(i,j) + dt * Phi(i,j) );
      % projection:
      Ptmp = max(0.0, Ptmp);
      Ptmp = min(Ptmp, Po(i,j));
      Pnew(i,j) = Ptmp;
    end
  end

  % W time step
  Wnew = W;  % copies unaltered initial b.c.s
  inputvol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      Wij = W(i,j);
      upe = up(alphV(i,j),  Wij,     W(i+1,j));
      upw = up(alphV(i-1,j),W(i-1,j),Wij);
      upn = up(betaV(i,j),   Wij,     W(i,j+1));
      ups = up(betaV(i,j-1), W(i,j-1),Wij);
      inputdepth = dt * Phi(i,j);
      dtlapW = mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
               muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1)));
      Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + dtlapW + ...
                  inputdepth;
      inputvol = inputvol + inputdepth * dA;
    end
  end

  % do volume accounting
  losevol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      if outline(i,j) < 0.5
        losevol = losevol + Wnew(i,j)*dA;
        Wnew(i,j) = 0.0;
      end
      % FIXME: if Phi<0 is possible, should check W<0 here
    end
  end

  % report on new state W
  sumnew = sum(sum(Wnew));
  volnew = sumnew * dA;
  vbalance = (volW + inputvol - losevol) - volnew;
  if volnew > 0, relbal = abs(vbalance/volnew); else relbal = nan; end
  fprintf('t = %.5f (a):  %7.3f  %7.3f        %.3f        %.0e\n',...
          (t+dt)/p.spera, max(max(Wnew)), sumnew/((Mx+1)*(My+1)),...
          volnew/(1e9),relbal)

  % actually update to new state W
  volW = volnew;
  W = Wnew;
  P = Pnew;
  t = t + dt;
end

   function z = up(v,L,R)
   if v >= 0
     z = v * L;
   else
     z = v * R;
   end

