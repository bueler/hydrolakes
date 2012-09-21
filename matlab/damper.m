function W = damper(x,y,b,h,magvb,floatmask,W0,Y0,Phi,ts,te,Nmin)
% DAMPER Runs a nearly-full-cavity Darcy flux subglacial hydrology
% model.
%
% ** DESCRIPTION **: This model combines an advection-diffusion
% equation for water thickness W with a penalized form of the
% "Y=W" full-cavity condition which gives evolution of the capacity
% thickness Y.
%
% ** INPUT DATA **:  The input data are given on (Mx+1) x (My+1) grids:
%   x = coordinate vector of length Mx+1
%   y = coordinate vector of length My+1
%   b = bed elevation; array of (Mx+1) x (My+1) values
%   h = surface elevation; <same size>
%   magvb = |v_b| = magnitude of sliding velocity; <same size>
%   floatmask = 1 if floating, 0 otherwise; <same size>
%   W0 = initial values for water thickness W; <same size>
%   Y0 = initial values for capacity thickness Y; <same size>
%   Phi = melt rate as m s-1; <same size>
% Note the ice thickness is  H = h - b  if  floatmask==0  and  H = 0
% otherwise.
%
% ** MODEL **:  The major ("state space") functions are W and Y.  Some
% derived functions are:
%   P_o = rhoi g H = overburden pressure
%   Z = tau^{-1} (Y-W) + c_1 (W_r - Y) |v_b|
%       /  0,                               if Z >= c_2 A Y P_o^3
%   P = |  P_o,                             if Z <= 0
%       \  P_o - Z^{1/3} (c_2 A Y)^{-1/3},  in all other cases
%   V = - (K / (rhow g)) grad P - K grad b
%   q = V W - K W grad W
% The hydraulic potential of the top of the water layer is
% psi = P + rhow g (b + W), but it turns out we do not need this internally.
% Note that the velocity has named components  V = (alpha,beta)  in the code.
% Also note the flux  q  is a sum of advective and diffusive components.
% The evolution equations are
%
%   dW/dt + div( V W ) = div ( K W grad W ) + Phi
%   dY/dt = - tau^{-1} (Y - W)
%
% Note that the equation for dY/dt is a damped (penalized) version of Y=W.
%
% ** USAGE **  The calling form is
%
%      W = damper(x,y,b,h,magvb,floatmask,W0,Y0,Phi,ts,te,Nmin)
%
% The input data has already been described.  Regarding the last three
% scalar arguments, the code runs from time  ts  to time  te  using at
% least  Nmin  time steps.  Thus the maximum time step is  (te-ts)/Nmin,
% but the code does adaptive time-stepping.  At each time step the code
% reports time step, time step restrictions (CFL for advection and
% criterion for diffusion), and total water amount.

% If floatmask==1 then W=0 at that point.  (FIXME: account for water lost.)

% the following parameters could be packaged into an input structure:
spera = 31556926.0;
rhoi = 910.0;   % kg m-3
rhow = 1028.0;  % kg m-3
r = rhoi / rhow;
g = 9.81;       % m s-2

% major model parameters:
K = 5.0e-3;            % m s-1   FIXME
tau = (1/50) * spera;  % s; = 1 week
Wr = 1.0;              % m
c1 = 0.500;            % m-1
c2 = 0.040;            % [pure]

c0 = K / (rhow * g);   % constant in velocity formula

% get grid parameters from W0; other fields must match but code does not check
[Mx, My] = size(W0);
Mx = Mx-1;  My = My-1;
dx = x(2) - x(1);
dy = y(2) - y(1);
% generate staggered grid
xs = 0.5 * (x(1:end-1) + x(2:end));
ys = 0.5 * (y(1:end-1) + y(2:end));

% bed slope components onto staggered grid
dbdx = (b(2:end,:) - b(1:end-1,:)) / dx;
dbdy = (b(:,2:end) - b(:,1:end-1)) / dy;

dtmax = (te - ts) / Nmin;

fprintf('running ...\n')
fprintf('                                ')
fprintf('        [max V (m/a)    dtCFL (a)    dtdiffusion (a)  -->  dt (a)]\n')
fprintf('t (a):   max W (m)  av W (m)  vol(10^6 km^3)  rel-bal\n')

t = ts;
W = W0;
Y = Y0;
volW = 0.0;
dA = dx*dy;

while t<te
  % pressure is nontrivial
  Po = rhoi * g * (h - b);
  Z = (1/tau) * (Y - W) + c1 * (Wr - Y) .* magvb;
  %     /  0,                               Z >= c_2 A Y P_o^3
  % P = |  P_o,                             Z <= 0
  %     \  P_o - Z^{1/3} (c_2 A Y)^{-1/3},  other cases
  P = Po;  % addresses Z <= 0 cases already, and allocates
  for i=1:Mx
    for j=1:My
      if Z(i,j) >= c2 * A * Y(i,j) * Po(i,j)^3
        P(i,j) = 0;
      elseif Z(i,j) > 0
        P(i,j) = Po(i,j) - Z(i,j)^(1/3) * (c2 * A * Y(i,j))^(-1/3);
      end
    end
  end

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

  nux = dt / dx;        nuy = dt / dy;
  mux = K * dt / dx^2;  muy = K * dt / dy^2;

  Wnew = W;  % copies unaltered zero b.c.s

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
      X = dt / tau;
      Ynew(i,j) = (1 / (1+x)) * Y(i,j) + (1 / (1 + (1/x))) * Wnew(i,j);
    end
  end

  losevol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      if floatmask(i,j) > 0.5
        losevol = losevol + Wnew(i,j)*dA;
        Wnew(i,j) = 0.0;
      end
    end
  end

  sumnew = sum(sum(Wnew));
  fprintf('t = %.5f (a):  %7.3f  %7.3f',...
          (t+dt)/spera, max(max(Wnew)), sumnew/((Mx+1)*(My+1)))

  volnew = sumnew * dA;
  balance = (volW + inputvol - losevol) - volnew;
  fprintf('        %.3f        %.0e\n',...
            volnew/(1e9*1e6),abs(balance/volnew))

  volW = volnew;
  W = Wnew;
  t = t + dt;
end

   function z = up(v,L,R)
   if v >= 0
     z = v * L;
   else
     z = v * R;
   end
