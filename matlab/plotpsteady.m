function plotpsteady;
% PLOTPSTEADY  Plot P = P(W) in steady state.

spera = 31556926.0;
rhoi  = 910.0;         % kg m-3
rhow  = 1028.0;        % kg m-3
g     = 9.81;          % m s-2

% major model parameters:
A  = 3.1689e-24;       % ice softness (Pa-3 s-1)
Wr = 1.0;              % m
c1 = 0.500;            % m-1
c2 = 0.040;            % [pure]
CC = c1 / (c2 * A)

h0 = 1000.0;         % m
Po = rhoi * g * h0;

W = 0.01:.01:1.1*Wr;

% four different amounts of sliding:
v0   = 100.0 / spera;  % m/s
Pno  = PofW(W,    0.0,Po,CC,Wr) / 1e5;
Plow = PofW(W, 0.1*v0,Po,CC,Wr) / 1e5;
Pmed = PofW(W, 1.0*v0,Po,CC,Wr) / 1e5;
Phgh = PofW(W,10.1*v0,Po,CC,Wr) / 1e5;

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

plot(W,Pno,W,Plow,W,Pmed,W,Phgh);
legend('no sliding; P=P_o  ','|v_b|=10 m/a','|v_b|=100 m/a','|v_b|=1000 m/a',...
       'location','southeast')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.1*Po/1e5])

  function P = PofW(W,vb,Po,CC,Wr)
  frac = CC * vb * max(0.0, Wr - W) ./ W;
  P = max(0.0, Po - frac.^(1/3));
  end

end

