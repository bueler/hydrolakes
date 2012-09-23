function W = conserve(x,y,b,h,floatmask,W0,Phi,ts,te,Nmin)
% CONSERVE Runs the minimal advection-diffusion model for subglacial
% hydrology from time  ts  to time  te  using at least  Nmin  time steps:
%   W_t + div q = Phi
%   q = - (K / rhow g) (grad psi) W
%   psi = P + rhow g (b + W)
%   P = s P_i = s rhoi g H
% where 
%   b = bed elevation
%   h = surface elevation
%   H = h - b  if  floatmask==0
% Note that
%   v = (alpha,beta) = - s K r grad h - K (1 - s r) grad b
% and the flux  q  is a sum of advective and diffusive components:
%   q = v W - K W grad W
% so the equation is
%   W_t + div (v W) = div (K W grad W) + Phi
% If floatmask==1 then W=0 at that point.  (FIXME: account for water lost.)
% Form:  W = conserve(x,y,b,h,floatmask,W0,Phi,ts,te,Nmin)

rhoi = 910.0;   % kg m-3
rhow = 1028.0;  % kg m-3
r = rhoi / rhow;
g = 9.81;       % m s-2
K = 1.0e-3;     % m s-1

spera = 31556926.0;

s = 1;  % s=1 is full overburden
sr = s * r;

[Mx, My] = size(W0);
dx = x(2) - x(1);
dy = y(2) - y(1);

xs = 0.5 * (x(1:end-1) + x(2:end));
ys = 0.5 * (y(1:end-1) + y(2:end));

dhdx  = (h(2:end,:)-h(1:end-1,:)) / dx;
dbdx  = (b(2:end,:)-b(1:end-1,:)) / dx;
alpha = - K * (sr * dhdx + (1-sr) * dbdx);
dhdy  = (h(:,2:end)-h(:,1:end-1)) / dy;
dbdy  = (b(:,2:end)-b(:,1:end-1)) / dy;
beta  = - K * (sr * dhdy + (1-sr) * dbdy);

%figure(5), scale=1.0e11;
%quiver(xs/1000,ys/1000,scale*alpha(:,1:end-1),scale*beta(1:end-1,:))
%title('velocity'), xlabel('x (km)'), ylabel('y (km)')

dtmax = (te - ts) / Nmin;

maxv = sqrt(max(max(alpha))^2 + max(max(beta))^2);
dtCFL = 0.5 / (max(max(alpha))/dx + max(max(beta))/dy);
fprintf('max |v| = %.3f m a-1  and  dtCFL = %.3f a\n',...
  maxv*spera, dtCFL/spera)

fprintf('running ...\n')
fprintf('                                ')
fprintf('max W (m)  av W (m)  vol(10^6 km^3)  rel-bal\n')

t = ts;
W = W0;
volW = 0.0;
dA = dx*dy;
while t<te
  dt = min(te-t,dtmax);

  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));

  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFF = 0.25 / (K * maxW * (1/dx^2 + 1/dy^2));
  dt = min([dt dtCFL dtDIFF]);
  fprintf('  t = %8.3f a [dt = %6.3f]:  ',(t+dt)/spera,dt/spera)

  nux = dt / dx;        nuy = dt / dy;
  mux = K * dt / dx^2;  muy = K * dt / dy^2;

  Wnew = W;  % copies unaltered zero b.c.s

  inputvol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
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
  fprintf('%7.3f  %7.3f',max(max(Wnew)),sumnew/(Mx*My))

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
