function [alpha,beta,W] = conserve(x,y,b,h,floatmask,W0,ts,te,Nmin)
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
%   q = v W - K W grad W.
% If floatmask==1 then W=0 at that point.  (FIXME: account for water lost.)

rhoi = 910.0;   % kg m-3
rhow = 1028.0;  % kg m-3
r = rhoi / rhow;
g = 9.81;       % m s-2
K = 1.0e-3;     % m s-1

s = 1;  % s=1 is full overburden

Phi = 0.01;   % 1 cm s-1

[Mx, My] = size(W0);



dtmax = (te - ts) / Nmin;
t = ts;
W = W0;

return  % FIXME

while t<te
  dt = min(te-t,dtmax);
  % Wea = east staggered, Wno = north staggered
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));
end

