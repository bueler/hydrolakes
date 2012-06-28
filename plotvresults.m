% PLOTVRESULTS  semilog plot to show how v depends on grid spacing

dxkm = [100 50 25 15 10 5];
vmaxpera = [694.514 2240.086 4378.39 6887.370 9992.122 12847.341];

semilogy(dxkm, vmaxpera, 'o', 'markersize', 14)
axis([0 105 500 20000])
xlabel('\Delta x  (km)')
ylabel('max |v|  (m a-1)')
set(gca,'xtick',[dxkm 0])
set(gca,'xticklabel',{'100','50','25','15','10','5','0'})
set(gca,'ytick',[1000 2000 5000 10000])
set(gca,'yticklabel',{'1000','2000','5000','10000'})
grid on

