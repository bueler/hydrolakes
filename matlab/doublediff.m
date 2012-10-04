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

FIXME FROM HERE

  % pressure is nontrivial
  Vcav = c1 * (Wr - Y) .* magvb;
  Vcav(Vcav < 0.0) = 0.0;
  if ward
    Z = (Y - W)/dt - (W - Wold)/dt + Vcav;
  else
    Z = (Y - W)/tau                + Vcav;
  end
  %     /  0,                               Z >= c_2 A Y P_o^3
  % P = |  P_o,                             Z <= 0
  %     \  P_o - Z^{1/3} (c_2 A Y)^{-1/3},  other cases
  P = Po;  % addresses Z <= 0 cases already, and allocates
  for i=1:Mx+1
    for j=1:My+1
      if i==1 | i==Mx+1 | j==1 | j==My+1  % boundary case
        P(i,j) = 0;
      elseif outline(i,j) < 0.5
        P(i,j) = 0;
      else                                % interior case
        gam = c2 * A * (Y(i,j) + Yeps);
        if Z(i,j) >= gam * Po(i,j)^3
          P(i,j) = 0;
        elseif Z(i,j) > 0
          P(i,j) = Po(i,j) - (Z(i,j) / gam)^(1/3);
        end
      end
    end
  end

  % FIXME:  bypass and go back to overburden
  %P = Po;

  % grad pressure
  dPdx = (P(2:end,:) - P(1:end-1,:)) / dx;
  dPdy = (P(:,2:end) - P(:,1:end-1)) / dy;

  % velocity  V = - c0 grad P - K grad b = (alpha,beta)
  alpha = - c0 * dPdx - K * dbdx;
  beta  = - c0 * dPdy - K * dbdy;
  %figure(5)
  %scale = 1.0e4 * spera;  % scale so that 10000 m/a
  %quiver(xs/1000,ys/1000,scale*alpha(:,1:end-1),scale*beta(1:end-1,:))
  %title('velocity'), xlabel('x (km)'), ylabel('y (km)')

  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));

  % determine time step adaptively
  maxv = sqrt(max(max(alpha))^2 + max(max(beta))^2);
  dtCFL = 0.5 / (max(max(alpha))/dx + max(max(beta))/dy);
  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFF = 0.25 / (K * maxW * (1/dx^2 + 1/dy^2));
  dt = min([te-t dtmax dtCFL dtDIFF]);
  fprintf('        [%.5e  %.5f  %.5f   -->  dt = %.5f (a)]\n',...
          maxv*spera, dtCFL/spera, dtDIFF/spera, dt/spera)

  X = dt / tau;  % appears in update formula for Y
  nux = dt / dx;        nuy = dt / dy;
  mux = K * dt / dx^2;  muy = K * dt / dy^2;

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
      Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + ...
          mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
          muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1))) + ...
          inputdepth;
      inputvol = inputvol + inputdepth * dA;

      % evolve capacity
      if ward
        %s0 = c1 * Wr * magvb(i,j);
        %s1 = c1 * magvb(i,j) + c2 * A * (Po(i,j) - P(i,j))^3;
        %Ynew(i,j) = ( Y(i,j) + s0 * dt ) / (1 + s1 * dt);
        Ynew(i,j) = 0.5 * ( Y(i,j) + 2 * Wnew(i,j) - W(i,j) );
      else
        Ynew(i,j) = (1 / (1+X)) * Y(i,j) + (1 / (1 + (1/X))) * Wnew(i,j);
      end
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

