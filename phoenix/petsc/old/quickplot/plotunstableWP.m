% PLOTUNSTABLEWP

close all
figure(1), clf
load alt.out
load kw.out
load por.out
plot(por(:,2),por(:,4),alt(:,2),alt(:,4),kw(:,2),kw(:,4))
xlabel('t (days)'), ylabel('max W (m)')
legend('porous medium equation','(W,P) system','(W,K) system')

