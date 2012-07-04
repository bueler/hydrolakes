function plotobstacle(N)
% PLOTOBSTACLE  Reference plot for obstacle problem.  See
% online note obstacleDOC.pdf.

if nargin<1, N=40; end

dx = 4.0 / N;
[xx,yy] = meshgrid(-2:dx:2,-2:dx:2);
rr = sqrt(xx.*xx+yy.*yy);

psi = - ones(size(xx));
psi(rr<=1) = sqrt(1 - rr(rr<=1) .* rr(rr<=1));

afree = 0.69797;
A = 0.68026;
B = 0.47152;
uu = zeros(size(xx));
uu(rr>afree) = -A * log(rr(rr>afree)) + B;
uu(rr<=afree) = psi(rr<=afree);

mesh(xx,yy,uu)
hold on, mesh(xx,yy,psi), hold off
xlabel x, ylabel y

