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
Pno  = psteady(p,Po,0.0,W) / 1e5;
% compare Flowers and Clarke function with W_crit = W_r
PFC  = Po * (W / p.Wr).^(7/2) / 1e5;

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

figure(1)
plot(W,Pno,'b',W,PFC,'k--')
hold on
sty = ['b'; 'g'; 'r'; 'c'];
aa = [0.1 1.0 10.0];
for j=1:3
  vb = aa(j) * v0;
  sbcube = p.c1 * vb / (p.c2 * p.A);
  Wc = sbcube * p.Wr - Po^3 * p.Y0;
  Wc = Wc / (sbcube + Po^3);
  Pc = psteady(p,Po,vb,Wc);
  fprintf('for vb = %8.3f m/a, Wc = %9.6f m and Pc = %10.3f Pa\n',...
          vb * p.spera, Wc, Pc)
  W = Wc:0.0005:1.2*p.Wr;
  P = psteady(p,Po,vb,W) / 1e5;
  plot(W,P,sty(j+1,:))
  plot(Wc,Pc,'k.','markersize',30)
end
hold off
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.3*Po/1e5])

figure(2)
bb = [2.0 1.0 0.5 0.2];  % H = 2000, 1000, 500, 200 m
for j=1:4
  myPo = bb(j) * Po;
  sbcube = p.c1 * v0 / (p.c2 * p.A);
  Wc = sbcube * p.Wr - myPo^3 * p.Y0;
  Wc = Wc / (sbcube + myPo^3);
  Pc = psteady(p,myPo,v0,Wc);
  W = Wc:0.0005:1.2*p.Wr;
  P = psteady(p,myPo,v0,W) / 1e5;
  plot(W,P,sty(j,:))
  hold on
  plot(Wc,Pc,'k.','markersize',30)
end
hold off
xlabel('W  (m)'), ylabel('P  (bar)')
axis([0 max(W) 0 1.15*2*Po/1e5])

end

