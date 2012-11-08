% REFINEVERIF call verifwater along a refinement path

M = [20 40 80 160 200];
dx = zeros(1,length(M));
err = ones(length(M),4);
for j = 1:length(M)
  tic
  dx(j) = 60e3 / M(j);
  err(j,:) = verifwater(0.1,M(j),0);
  toc
end

loglog(dx,err(1,:),'*')
grid on
xlabel('\Delta x  (m)')
