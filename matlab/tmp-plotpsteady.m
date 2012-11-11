function plotpsteady;
% PLOTPSTEADY  Plot P = P(W) in steady state.  Generates figures for dampnotes.pdf.

p = params();

CC = p.c1 / (p.c2 * p.A);

h0 = 1000.0;         % m
Po = p.rhoi * p.g * h0;

W = 0.0:0.0005:1.2*p.Wr;

% four different amounts of sliding:
v0   = 100.0 / p.spera;  % m/s
Pno  = PofW(p,W,    0.0,Po) / 1e5;
Plow = PofW(p,W, 0.1*v0,Po) / 1e5;
Pmed = PofW(p,W, 1.0*v0,Po) / 1e5;
Phgh = PofW(p,W,10.0*v0,Po) / 1e5;

% compare Flowers and Clarke function
WcritFC = 0.5;
PFC  = Po * (W / WcritFC).^(7/2);
PFC  = PFC / 1e5;

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

figure(1)
plot(W,Pno,W,Plow,W,Pmed,W,Phgh,W,PFC,'k--');
legend('no sliding; P=P_o  ','|v_b|=10 m/a','|v_b|=100 m/a','|v_b|=1000 m/a',...
       'location','southeast')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.3*Po/1e5])
text(0.05,108,'(a)','fontsize',18)

figure(2)
PP1 = PofW(p,W,v0,0.2*Po) / 1e5; % H =  200 m
PP2 = PofW(p,W,v0,0.5*Po) / 1e5; % H =  500 m
PP3 = Pmed;                      % H = 1000 m
PP4 = PofW(p,W,v0,2*Po) / 1e5;   % H = 2000 m
plot(W,PP1,W,PP2,W,PP3,W,PP4);
legend('H = 200 m','H = 500 m','H = 1000 m','H = 2000 m',...
       'location','Best')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.15*2*Po/1e5])
text(0.05,1.08*2*Po/1e5,'(b)','fontsize',18)

  function P = PofW(p,W,vb,Po)
    CC = p.c1 / (p.c2 * p.A);
    frac = CC * vb * max(0.0, p.Wr - W) ./ (W + p.W0);
    P = Po + p.N0 - frac.^(1/3);
    P(P<0.0) = 0.0;
    P(P>Po) = Po;
  end

end

