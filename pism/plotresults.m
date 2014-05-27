% PLOTRESULTS  read verifTestP.results and plot

load('verifTestP.results')
dx = 50.0e3 ./ (verifTestP(:,1) - 1);  % in m
averrW  = verifTestP(:,2);
maxerrW = verifTestP(:,3);
averrP  = verifTestP(:,4) / 1.0e5;  % in bar
maxerrP = verifTestP(:,5) / 1.0e5;  % in bar

format short g
pW = polyfit(log(dx(1:end-1)),log(averrW(1:end-1)),1);
pP = polyfit(log(dx(1:end-1)),log(averrP(1:end-1)),1);
fprintf('fit to all average error data except finest grid\n')
fprintf('  (av err bwat) = O(dx^{%.2f}),  (av err bwp) = O(dx^{%.2f})\n',pW(1),pP(1))

figure(1)
subplot(4,1,1:2)
loglog(dx,averrW,'ks',dx,maxerrW,'kx','markersize',10)
grid on,  ylabel('error in W  (m)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(averrW) 1.1*max(maxerrW)])
set(gca,'XTick',flipud(dx))
set(gca,'XTickLabel',{'','','','',''})
set(gca,'YTick',[0.001 0.01])
set(gca,'YTickLabel',{'0.001','0.01'})
legend('average W error','maximum W error','Location','NorthWest')

subplot(4,1,3:4)
loglog(dx,averrP,'ks',dx,maxerrP,'kx','markersize',10)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in P  (bar)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(averrP) 1.1*max(maxerrP)])
set(gca,'XTick',flipud(dx))
set(gca,'XTickLabel',{'125','250','500','1000','2000'})
set(gca,'YTick',[0.01 0.1])
set(gca,'YTickLabel',{'0.01','0.1'})
legend('average P error','maximum P error','Location','NorthWest')

set(gcf,'Position',[100,100,600,600])

print -dpdf refineWPpism.pdf

