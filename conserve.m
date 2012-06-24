function [W,xx,yy,qx,qy] = conserve(x,y,b,h,floatmask,W0,ts,te,Nmin)
% CONSERVE Runs the minimal advection-diffusion model for subglacial
% hydrology from time  ts  to time  te  using at least  Nmin  time
% steps:
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
% Form:  [W,xx,yy,qx,qy] = conserve(x,y,b,h,floatmask,W0,ts,te,Nmin)
% Outputs  xx,yy,qx,qy  are only returned if desired.  They are suitable
% for use with quiver (

rhoi = 910.0;   % kg m-3
rhow = 1028.0;  % kg m-3
r = rhoi / rhow;
g = 9.81;       % m s-2
K = 1.0e-3;     % m s-1

spera = 31556926.0;

s = 1;  % s=1 is full overburden
sr = s * r;

Phi = 0.01 / spera;   % 1 cm a-1

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
fprintf('running ...\n')
fprintf('  [max |v| = %.3f m a-1  and  dtCFL = %.3f a]\n',...
  maxv*spera, dtCFL/spera)

fprintf('                                       ')
fprintf('    max W     av W\n')

t = ts;
W = W0;
while t<te
  dt = min(te-t,dtmax);

  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));

  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFF = 0.25 / (K * maxW * (1/dx^2 + 1/dy^2));
  dt = min([dt dtCFL dtDIFF]);
  fprintf('  t = %8.3f a  and  dt = %6.3f a:  ',t/spera,dt/spera)

  nux = dt / dx;        nuy = dt / dy;
  mux = K * dt / dx^2;  muy = K * dt / dy^2;

  Wnew = W;  % copies unaltered zero b.c.s
  for i=2:length(x)-1
    for j=2:length(y)-1
      if floatmask(i,j)<0.5
        Wij = W(i,j);
        upe = up(alpha(i,j),  Wij,     W(i+1,j));
        upw = up(alpha(i-1,j),W(i-1,j),Wij);
        upn = up(beta(i,j),   Wij,     W(i,j+1));
        ups = up(beta(i,j-1), W(i,j-1),Wij);
        Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + ...
            mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
            muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1))) + ...
            dt * Phi;
      else
        Wnew(i,j) = 0.0;
      end
    end
  end
  W = Wnew;
  fprintf(' %8.3f  %8.3f\n',max(max(W)),sum(sum(W))/(Mx*My))
  t = t + dt;
end

if nargout>1
  xx = 0.5*(xs(1:end-1) + xs(2:end));
  yy = 0.5*(ys(1:end-1) + ys(2:end));
  qx = zeros(Mx-2,My-2); %FIXME
  qy = zeros(Mx-2,My-2); %FIXME
end

   function z = up(v,L,R)
   if v >= 0
     z = v * L;
   else
     z = v * R;
   end
