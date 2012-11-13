function plotpsteady;
% PLOTPSTEADY  Plot P = P(W) in steady state with various assumptions
% about the overburden pressure Po and the sliding speed |v_b|.
% Generates figures for dampnotes.pdf.

p = params();
CC = p.c1 / (p.c2 * p.A);

h0 = 1000.0;             % m
Po = p.rhoi * p.g * h0;  % Pa; base overburden
v0   = 100.0 / p.spera;  % m/s; base sliding rate

% full range of W for vb=0 and F&C cases
W = 0.0:0.0005:1.2*p.Wr;
Pno  = PofW(p,W,    0.0,Po) / 1e5;
% compare Flowers and Clarke function with W_crit = W_r
PFC  = Po * (W / p.Wr).^(7/2) / 1e5;

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

figure(1)
plot(W,Pno,'k:',W,PFC,'k--')
hold on
sty = ['ks-'; 'kd-'; 'k*-']
aa = [0.1 1.0 10.0]
for j=1:3
  vb = aa(j) * v0;
  sbcube = p.c1 * vb / (p.c2 * p.A);
  Wc = sbcube * p.Wr - Po^3 * p.Y0;
  Wc = Wc / (sbcube + Po^3);
  Pc = PofW(p,Wc,vb,Po);
  plot(Wc,Pc,'k.','markersize',18)
  fprintf('for vb = %8.3f m/a, Wc = %9.6f m and Pc = %10.3f Pa\n',
          vb * p.spera, Wc, Pc)
  W = Wc:0.0005:1.2*p.Wr;
  P = PofW(p,W,vb,Po) / 1e5;
  plot(W,P,sty(j,:))
end
hold off
legend('no sliding; P=P_o  ','|v_b|=10 m/a','|v_b|=100 m/a','|v_b|=1000 m/a',...
       'location','southeast')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.3*Po/1e5])
text(0.05,108,'(a)','fontsize',18)

figure(2)
W = 0.0:0.0005:1.2*p.Wr;
PP1 = PofW(p,W,v0,0.2*Po) / 1e5; % H =  200 m
PP2 = PofW(p,W,v0,0.5*Po) / 1e5; % H =  500 m
PP3 = Pmed;                      % H = 1000 m
PP4 = PofW(p,W,v0,2*Po) / 1e5;   % H = 2000 m
plot(W,PP1,W,PP2,W,PP3,W,PP4);
legend('H = 200 m','H = 500 m','H = 1000 m','H = 2000 m',...
       'location','east')
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.15*2*Po/1e5])
text(0.05,1.08*2*Po/1e5,'(b)','fontsize',18)

  function P = PofW(p,W,vb,Po)
    sbcube = p.c1 * vb / (p.c2 * p.A);
    frac = max(0.0, p.Wr - W) ./ (W + p.Y0);
    P = Po - (sbcube * frac).^(1/3);
    P(P<0.0) = 0.0;
    P(P>Po) = Po;
  end

end

