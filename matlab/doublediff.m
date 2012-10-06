function [W, P] = doublediff(x,y,b,h,magvb,outline,W0,P0,Phi,ts,te,Nmin)
% DOUBLEDIFF  A nearly-full-cavity, Darcy flux subglacial hydrology model.
% This model combines an advection-diffusion equation for water thickness W
% with a penalized form of the "Y=W" full-cavity condition.  The latter
% evolves the hydraulic potential psi by a diffusion equation.
%
% ** INPUT DATA **:  The input data are given on (Mx+1) x (My+1) grids:
%   x = coordinate vector of length Mx+1
%   y = coordinate vector of length My+1
%   b = bed elevation; array of (Mx+1) x (My+1) values
%   h = surface elevation; <same size>
%   magvb = |v_b| = magnitude of sliding velocity; <same size>
%   outline = 1 if inside modeled region, 0 otherwise; <same size>
%   W0 = initial values for water thickness W; <same size>
%   P0 = initial values for pressure P; <same size>
%   Phi = melt rate as m s-1; <same size>
% Note that where outline==0 we set the water thickness W to zero.
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
% The input data has already been described.  Regarding the last three
% scalar arguments, the code runs from times  ts  to  te  using at
% least  Nmin  time steps.  Thus the maximum time step is  (te-ts)/Nmin,
% but the code does adaptive time-stepping.  At each time step the code
% reports time step, time step restrictions (CFL for advection and
% criterion for diffusion), and total water amount.  The output variables
% are the final values of the state variables W,P.

spera = 31556926.0;
rhoi = 910.0;   % kg m-3
rhow = 1028.0;  % kg m-3
r = rhoi / rhow;
g = 9.81;       % m s-2

% major model parameters:
A = 3.1689e-24;        % ice softness (Pa-3 s-1)
K = 1.0e-2;            % m s-1   FIXME: want Kmax or Kmin according to W > Wr
tau = (1/50) * spera;  % s; about 1 week
Wr = 1.0;              % m
c1 = 0.500;            % m-1
c2 = 0.040;            % [pure]

E0 = 1.0;         % m  FIXME:  not small; what is right

c0 = K / (rhow * g);   % constant in velocity formula

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

dtmax = (te - ts) / Nmin;

t = ts;
W = W0;
P = P0;
volW = 0.0;
dA = dx*dy;

fprintf('running ...\n\n')
fprintf('        [max V (m/a)    dtCFL (a)    dtDIFF_W (a)   dtDIFF_P (a)  -->  dt (a)]\n')
fprintf('t (a):   max W (m)  av W (m)  vol(10^6 km^3)  rel-bal\n\n')

while t<te
  % overburden pressure
  Po = rhoi * g * (h - b);
  Po(h <= 0.0) = 0.0;   % if no ice, ignor depth to bed

  % terms in pressure equation
  Ocav = c1 * (Wr - W) .* magvb;
  Ocav(Ocav < 0.0) = 0.0;
  Ccrp = c2 * (Po - P).^3 .* W;

FIXME

  % grad pressure and grad water
  dPdx = (P(2:end,:) - P(1:end-1,:)) / dx;
  dPdy = (P(:,2:end) - P(:,1:end-1)) / dy;

  % velocity  V = - c0 grad P - K grad b = (alpha,beta)
  alpha = - c0 * dPdx - K * dbdx;
  beta  = - c0 * dPdy - K * dbdy;

  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));

  % determine time step adaptively
  maxv = sqrt(max(max(alpha))^2 + max(max(beta))^2);
  dtCFL = 0.5 / (max(max(alpha))/dx + max(max(beta))/dy);
  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFFW = 0.25 / (K * maxW * (1/dx^2 + 1/dy^2));
  maxPo = max(max(Po));
  dtDIFFP = 0.25 / ((c0 * maxW * maxPo / E0) * (1/dx^2 + 1/dy^2));
  dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  fprintf('        [%.5e  %.5f  %.5f  %.5f   -->  dt = %.5f (a)]\n',...
          maxv*spera, dtCFL/spera, dtDIFFW/spera, dtDIFFP/spera, dt/spera)

  X = dt / tau;  % appears in update formula for Y
  nux = dt / dx;
  nuy = dt / dy;
  mux = K * dt / dx^2;
  muy = K * dt / dy^2;
  pux = (c0 / E0) * dt / dx^2;
  puy = (c0 / E0) * dt / dy^2;

FIXME from here

  Wnew = W;  % copies unaltered initial b.c.s
  Ynew = Y;  % copies unaltered initial b.c.s

  inputvol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      % conserve water
      Wij = W(i,j);
      upe = up(alpha(i,j),  Wij,     W(i+1,j));
      upw = up(alpha(i-1,j),W(i-1,j),Wij);
      upn = up(beta(i,j),   Wij,     W(i,j+1));
      ups = up(beta(i,j-1), W(i,j-1),Wij);
      inputdepth = dt * Phi(i,j);
      dtlapW = mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
               muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1)));
      Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + lapW + ...
          inputdepth;
      inputvol = inputvol + inputdepth * dA;

FIXME
    end
  end

  losevol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      if outline(i,j) < 0.5
        losevol = losevol + Wnew(i,j)*dA;
        Wnew(i,j) = 0.0;
      end
    end
  end

  % report on new state W,Y
  sumnew = sum(sum(Wnew));
  volnew = sumnew * dA;
  balance = (volW + inputvol - losevol) - volnew;
  fprintf('t = %.5f (a):  %7.3f  %7.3f        %.3f        %.0e\n',...
          (t+dt)/spera, max(max(Wnew)), sumnew/((Mx+1)*(My+1)),...
          volnew/(1e9*1e6),abs(balance/volnew))

  if ward
    Wold = W;
  end

  % actually update to new state W,Y
  volW = volnew;
  W = Wnew;
  Y = Ynew;
  t = t + dt;
end

   function z = up(v,L,R)
   if v >= 0
     z = v * L;
   else
     z = v * R;
   end

