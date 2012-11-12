function p = params(N0,W0,E0);

p.spera = 31556926.0;
p.rhoi  = 910.0;         % kg m-3
p.rhow  = 1000.0;        % kg m-3
p.g     = 9.81;          % m s-2

% major model parameters:
p.A  = 3.1689e-24;       % ice softness (Pa-3 s-1)
p.K  = 1.0e-2;           % m s-1   FIXME: want Kmax or Kmin according to W > Wr
p.Wr = 1.0;              % m
p.c1 = 0.500;            % m-1
p.c2 = 0.040;            % [pure]
p.c0 = p.K / (p.rhow * p.g);   % constant in velocity formula

% regularizations of closing term
if nargin < 1
  p.N0 = p.rhoi * p.g * 1.0; % Pa
else
  p.N0 = N0;
end
if nargin < 2
  p.W0 = 0.01;            % m
else
  p.W0 = W0;
end

% regularization of "W=Y" pressure closure
if nargin < 1
  p.E0 = 1.0;            % m; what is optimal?
else
  p.E0 = E0;
end

