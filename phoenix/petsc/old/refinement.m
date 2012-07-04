function refinement
% REFINEMENT  Produces a figure with
%   dx  versus  L^1, L^\infty error
% for runs in which dt scales linearly with dx

set(0,'defaultlinelinewidth',2.0)

L = 40000;
Wnorm1 = 5.65232e+07;  % m^3; from 201x201 calculation (good enough for normalization)

C=load('refine_1517.txt')
C = reshape(C,4,7)';
dtdays = 20 ./ C(:,1);
dx = L ./ (C(:,2)-1);
p1 = polyfit(log(dx(3:7)),log(C(3:7,3)/Wnorm1),1)
pinf = polyfit(log(dx(3:7)),log(C(3:7,4)),1)
loglog(dx,C(:,3)/Wnorm1,'bo','markersize',12,...
       dx,C(:,4),'r*','markersize',12)
xlabel('dx  (m)'), grid on
set(gca,'xtick',dx)
set(gca,'xticklabel',{'2000','1000','500','250','125','63','31'})
set(gca,'ytick',10.^(-3.5:.5:-1.5))
set(gca,'yticklabel',{'10^{-3.5}','10^{-3}','10^{-2.5}','10^{-2}','10^{-1.5}'})
%legend('|W-W_{exact}|_1 / |W|_1','|W-W_{exact}|_{inf}','location','southeast')
hold on
loglog(dx,exp(p1(1)*log(dx)+p1(2)),'b:','linewidth',1.0)
loglog(dx,exp(pinf(1)*log(dx)+pinf(2)),'r:','linewidth',1.0)
hold off
axis([31.25 2000 10^(-3.7) 10^(-1.3)])
end

