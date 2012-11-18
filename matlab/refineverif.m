% REFINEVERIF call verifwater along a refinement path

tyears = 1/12;  % one month coupled runs

%M = [10 20 40 80 160 200];
M = [10 20 40 80 160];
%M = [10 20 40];
dx = zeros(1,length(M));
err = ones(length(M),4);
for j = 1:length(M)
  tic
  dx(j) = 60e3 / M(j);
  err(j,:) = verifwater(tyears,M(j),false);
  toc
end

figure(21)
loglog(dx,err(:,1),'ko',dx,err(:,3),'k*','markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in W  (m)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(err(:,1)) 1.1*max(err(:,3))])
set(gca,'XTick',fliplr(dx))
legend('average W error','maximum W error')

% put P errors in bar
err(:,2) = err(:,2)/1e5;
err(:,4) = err(:,4)/1e5;
figure(22)
loglog(dx,err(:,2),'ko',dx,err(:,4),'k*','markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in P  (bar)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(err(:,2)) 1.1*max(err(:,4))])
set(gca,'XTick',fliplr(dx))
legend('average P error','maximum P error')

pW = polyfit(log(dx(2:end)),log(err(2:end,1)'),1);
% pW =   1.0534  -12.4881
pP = polyfit(log(dx(2:end)),log(err(2:end,2)'),1);
% pP =   0.9556   -9.4883
%FIXME: add to figures?

