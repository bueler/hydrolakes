% REFINEVERIF call verifwater along a refinement path

tyears = 1/12;  % one month coupled runs
reload = false;

M = [10 20 40 80 160];
%M = [10 20 40 80];
%M = [10 20 40];

if reload
    load('PWerrdx.mat')
    whos
else
    dx = zeros(1,length(M));
    err = ones(2*length(M),4);
    for j = 1:length(M)
      dx(j) = 60e3 / M(j);
      tic
      err(2*j-1,:) = verifwater(tyears,M(j),false,true);  % limiter
      toc, tic
      err(2*j,:)   = verifwater(tyears,M(j),false,false); % 1st-order
      toc
    end
    save('PWerrdx.mat')
end

figure(21)
loglog(dx,err(1:2:2*length(M)-1,1),'ko-',dx,err(2:2:2*length(M),1),'ks',...
       dx,err(1:2:2*length(M)-1,3),'k*-',dx,err(2:2:2*length(M),3),'kx',...
       'markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in W  (m)')
axis([0.9*min(dx) 1.1*max(dx) 0.9*min(err(:,1)) 1.1*max(err(:,3))])
set(gca,'XTick',fliplr(dx))
legend('average W error (limiter)','average W error (first-order)',...
       'maximum W error (limiter)','maximum W error (first-order)',...
       'Location','SouthEast')
print -dpdf refineW.pdf

% put P errors in bar
err(:,2) = err(:,2)/1e5;
err(:,4) = err(:,4)/1e5;
figure(22)
loglog(dx,err(1:2:2*length(M)-1,2),'ko-',dx,err(2:2:2*length(M),2),'ks',...
       dx,err(1:2:2*length(M)-1,4),'k*-',dx,err(2:2:2*length(M),4),'kx',...
       'markersize',12)
grid on,  xlabel('\Delta x  (m)'), ylabel('error in P  (bar)')
%axis([0.9*min(dx) 1.1*max(dx) 0.9*min(err(:,2)) 1.1*max(err(:,4))])
axis([0.9*375 1.1*6000 0.009 1.05])
set(gca,'XTick',fliplr(dx))
legend('average P error (limiter)','average P error (first-order)',...
       'maximum P error (limiter)','maximum P error (first-order)',...
       'Location','SouthEast')
print -dpdf refineP.pdf

pWlimiter = polyfit(log(dx(2:end)),log(err(3:2:2*length(M)-1,1)'),1)
pWfirst   = polyfit(log(dx(2:end)),log(err(4:2:2*length(M),  1)'),1)
pPlimiter = polyfit(log(dx(2:end)),log(err(3:2:2*length(M)-1,2)'),1)
pPfirst   = polyfit(log(dx(2:end)),log(err(4:2:2*length(M),  2)'),1)

