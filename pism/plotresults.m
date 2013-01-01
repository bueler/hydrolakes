% PLOTRESULTS  read verifTestP.results and plot

load('verifTestP.results')
dx = 50.0e3 ./ (verifTestP(:,1) - 1);  % in m
averrW  = verifTestP(:,2);
maxerrW = verifTestP(:,3);
averrP  = verifTestP(:,4) / 1.0e5;  % in bar
maxerrP = verifTestP(:,5) / 1.0e5;  % in bar

figure(1)
loglog(dx,averrW,'ks',dx,maxerrW,'kx','markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in W  (m)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(averrW) 1.1*max(maxerrW)])
set(gca,'XTick',flipud(dx))
legend('average W error','maximum W error','Location','SouthEast')
print -dpdf refineWpism.pdf

figure(2)
loglog(dx,averrP,'ks',dx,maxerrP,'kx','markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in P  (bar)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(averrP) 1.1*max(maxerrP)])
set(gca,'XTick',flipud(dx))
legend('average P error','maximum P error','Location','SouthEast')
print -dpdf refinePpism.pdf

format short g
pW = polyfit(log(dx(2:end)),log(averrW(2:end)),1)
pP = polyfit(log(dx(2:end)),log(averrP(2:end)),1)

