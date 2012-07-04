function newtonconverge
% NEWTONCONVERGE   Run newtonconverge.sh before this one and generate
% ascii file converge1.txt.  This function produces a figure
%   dt  versus  residual norm 
% for one time step.

set(0,'defaultlinelinewidth',2.0)

C = load('converge1.txt');
i21  = find(C(:,1)==21);
i81  = find(C(:,1)==81);
i321 = find(C(:,1)==321);
loglog(C(i21,3),C(i21,4),'bo','markersize',12)
grid on, hold on
dtc = getcritical( 20)
loglog([dtc dtc], [1e-20 1e-4],'b-')
dte = getexplicit( 20)
loglog([dte dte], [1e-20 1e-4],'b--','linewidth',3.0)
loglog(C(i81,3),C(i81,4),'r*','markersize',12)
dtc = getcritical( 80)
loglog([dtc dtc], [1e-20 1e-4],'r-')
dte = getexplicit( 80)
loglog([dte dte], [1e-20 1e-4],'r--','linewidth',3.0)
loglog(C(i321,3),C(i321,4),'g+','markersize',12)
dtc = getcritical(320)
loglog([dtc dtc], [1e-20 1e-4],'g-')
dte = getexplicit(320)
loglog([dte dte], [1e-20 1e-4],'g--','linewidth',3.0)
hold off
xlabel('dt  (days)'), ylabel('2-norm of residual')
axis([10^(-2.5) 10^2.5 1e-20 1e-4])

  function dtc = getcritical(M)
  L = 40000;  % m
  maxD = 35.731;  % m2s-1; from max D for Mx=My=21 case at tstart=10 days
  secperday = 3600.0 * 24.0;
  dx = L / M;
  dtc = 3.5 * L^0.9 * dx^1.1 / maxD;  % speculative, but right units
  dtc = dtc / secperday;
  end

  function dtexpl = getexplicit(M)
  L = 40000;  % m
  maxD = 35.731;  % m2s-1; from max D for Mx=My=21 case at tstart=10 days
  secperday = 3600.0 * 24.0;
  dx = L / M;
  dtexpl = dx^2 / maxD;
  dtexpl = dtexpl / secperday;
  end

end
