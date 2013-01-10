function Kcompare;
% Kcompare  Plot K(W).

p = params();
K0 = 0.001;
K1 = 0.01;

W = 0.0:0.01:5*p.Wr;

K = K1 + (K0-K1) * exp(-max(0*W,W-p.Wr)/p.Wr);

set(0,'defaultlinelinewidth',3.0)
set(0,'defaultaxesfontsize',16.0)

figure(1)
plot(W,K,'k',W,K0*ones(size(W)),'k--',W,K1*ones(size(W)),'k:')
grid on
hold on
text(3.0,0.0018,'K(W)=K_0','FontSize',18.0)
text(2.2,0.0065,'K(W)','FontSize',18.0)
text(0.7,0.0104,'K_1','FontSize',18.0)
hold off
axis([min(W) max(W) 0.0 1.2*K1])
xlabel('W  (m)')
ylabel('hydraulic conductivity  (m s-1)')

