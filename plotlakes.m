% PLOTLAKES  load and minimal plot of lakes from Siegert (2005)

A = load('siegert_2005_lakes.txt');
B = load('vostok.txt');

plot(A(:,1)/1000,A(:,2)/1000,'o','markersize',14)
hold on
plot(B(:,2)/1000,B(:,3)/1000,'r.')
hold off
xlabel('x  (km)'), ylabel('y  (km)')
axis equal
