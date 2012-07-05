% PLOTLAKES  load and minimal plot of lakes from Siegert (2005)

A = load('siegert_2005_lakes.txt');
B = load('vostok.txt');

plot(A(:,1)/1000,A(:,2)/1000,'ro','markersize',7,'linewidth',2)
hold on
plot(B(:,2)/1000,B(:,3)/1000,'g.','markersize',7)
hold off
%xlabel('x  (km)'), ylabel('y  (km)')
axis([-2800.000 3195.000 -2800.000 3195.000])
axis equal
set(gca,'xtick',[])
set(gca,'ytick',[])

