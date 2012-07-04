function [xx,yy,W] = baren(frames,delay)
% BAREN  Plot the Barenblatt solution of the vertically-
% integrated porous medium equation.  See draft PISM hydrology
% model based on Flowers & Clarke (2002).
% Example 1: default movie
%   >> baren;
% Example 2: show first five frames with big delay
%   >> baren(5,1.0);
% Example 3: contour plot of fifth frame
%   >> [xx,yy,W] = baren(5,0.0);
%   >> [c,h] = contour(xx/1000,yy/1000,W*100,1:2:15);
%   >> axis equal, clabel(c,h);
%   >> title('W in cm'), xlabel('x (km)'), ylabel('y (km)')

if nargin < 1, frames = 20; end
if nargin < 2, delay = 0.1; end

% values from Flowers & Clarke (2002):
p.sigma = 7/2;
p.rhoi  = 910.0;  % kg m-3; ice density; also PISM default
p.rhow  = 1000.0; % kg m-3; freshwater density; also PISM default
p.Kmax  = 0.01;   % m s-1;  maximum of hydraulic conductivity
p.Kmin  = 0.0001; % m s-1;  minimum of hydraulic conductivity

% constant hydraulic conductivity in q=0 case:
K  = sqrt(p.Kmax*p.Kmin)  % geometric average of Kmax & Kmin

% Barenblatt solution parameters
Y0 = 1.0        % m;      critical water thickness
H0 = 3000       % m;      typical ice sheet thickness
W0 = 200        % m;      center water thickness in Barenblatt   
t0 = ((p.sigma+1) * p.rhow * Y0^p.sigma) ... %  Baren. timescale
       / (p.sigma * K * p.rhoi * H0 * W0^(p.sigma-2))
R0 = 2 * (p.sigma+1) * W0 / sqrt(p.sigma)  % wetting radius at t0
V0 = 4 * pi * (p.sigma+1) * W0^3

% setup time axis, and inform user
tstart = t0 * 10^6;
tend = t0 * 10^11;
if frames == 1, tend = tstart; end
tlist  = linspace(tstart,tend,frames);
secperday = 3600 * 24;
fprintf('BAREN: movie of Barenblatt solution as applied to\n')
fprintf('  subglacial hydrology model; shows %d frames\n',...
        frames)
fprintf('  from    %.0e t0 = %12.5f days\n  to      %.0e t0 = %12.5f days (=%.3f years)\n',...
        tstart/t0, tstart/secperday, tend/t0, tend/secperday, ...
        tend/(365.24*secperday))
if frames > 1
  fprintf('  with %.3f day time-steps between frames\n',...
        (tlist(2)-tstart)/secperday)
end

% grid on which to display
L = 2.0e4;   % 20 km so region is 40km x 40km
dx = 200;    % 200 m
x = -L:dx:L;  [xx, yy] = meshgrid(x);
rsqr = xx.*xx + yy.*yy;

% setup color axis
W0start = W0 * (t0/tstart)^(2/9);
W0end   = W0 * (t0/tend)^(2/9);
cmap = colormap;  cmap(1,:) = [1 1 1];  % show W=0 as white
colormap(cmap)
fprintf('  color shows  log_{10} W(t,r)  with W in meters\n')

% compute Barenblatt solution and show movie
for m = 1:frames
  t   = tlist(m);
  trp = (t0 / t)^(1/(p.sigma+1));
  A   = (p.sigma * trp) / (4 * (p.sigma+1)^2 * W0^2);
  psi = 1 - A * rsqr;
  psi = max(psi,0.0);
  W   = W0 * trp * psi.^(1/p.sigma);
  imagesc(x/1000,x/1000,log10(W),...
          [log10(0.0001*W0end) log10(W0start)])
  colorbar,  axis equal
  ylabel('y  (km)'),  xlabel('x  (km)')
  title(sprintf('W at t=%.6f  (days);   W(t,0) = %.6f  (m)',...
                t / secperday, W0 * trp))
  % optional:  print(sprintf('frame%03d.pdf',m),'-dpdf') 
  pause(delay)
end

