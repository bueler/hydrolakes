function plotpsteady;
% PLOTPSTEADY  Plot P = P(W) in steady state.

p = params();

CC = p.c1 / (p.c2 * p.A);

h0 = 1000.0;         % m
Po = p.rhoi * p.g * h0;

W = 0.0:0.0005:1.2*p.Wr;

% four different amounts of sliding:
v0   = 100.0 / p.spera;  % m/s
Pno  = PofW(W,    0.0,Po,CC,p.Wr) / 1e5;
Plow = PofW(W, 0.1*v0,Po,CC,p.Wr) / 1e5;
Pmed = PofW(W, 1.0*v0,Po,CC,p.Wr) / 1e5;
Phgh = PofW(W,10.1*v0,Po,CC,p.Wr) / 1e5;

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

plot(W,Pno,W,Plow,W,Pmed,W,Phgh);
legend('no sliding; P=P_o  ','|v_b|=10 m/a','|v_b|=100 m/a','|v_b|=1000 m/a',...
       'location','southeast')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.1*Po/1e5])

  function P = PofW(W,vb,Po,CC,Wr)
    if CC*vb > 0
      frac = CC * vb * max(0.0, Wr - W(W>0)) ./ W(W>0);
      P(W>0) = max(0.0, Po - frac.^(1/3));
      P(W<=0)= 0;
    else
      P = Po * ones(size(W));
    end
  end

end

