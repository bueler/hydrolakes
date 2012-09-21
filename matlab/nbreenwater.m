% NBREENWATER  **FIXME this comment**
% Calls:   BUILDNBREEN, CONSERVEWATER, DAMPER
% Depends: NETCDF

% % if nargin<1, tyears=100.0; end
% % if nargin<2, nn=10; end

tyears=5.0;
nn=4;

filename = 'nbreen_input.nc';

fprintf('reading variables x,y,thk,topg,usurf,icemask\n  from text files %s\n',filename)
[x,y,thk,topg,usurf,icemask, outline] = buildnbreen(0,filename);

fprintf('x range: %.3f km to %.3fkm\n',min(x)/1000,max(x)/1000)
fprintf('y range: %.3f km to %.3fkm\n',min(y)/1000,max(y)/1000)

dx = x(2)-x(1);
fprintf('subsampling to %d m resolution\n', dx*nn)
clear dx
x = x(1:nn:end);
y = y(1:nn:end);
thk = thk(1:nn:end,1:nn:end);
topg = topg(1:nn:end,1:nn:end);
usurf = usurf(1:nn:end,1:nn:end);
icemask = icemask(1:nn:end,1:nn:end);
outline = outline(1:nn:end,1:nn:end);

fprintf('showing topg, usurf, icemask\n')

figure(1)
imagesc(x/1000,flipud(-y)/1000,flipud(topg)), colorbar
title('bed elevation  (m)')
xlabel('x (km)'), ylabel('y (km)')
%set(gcf,'PaperPositionMode','auto')
%print('-dpsc2',strcat('fig_bed.eps'))

figure(2)
imagesc(x/1000,flipud(-y)/1000,flipud(usurf)), colorbar
title('surface elevation  (m)')
xlabel('x (km)'), ylabel('y (km)')
%set(gcf,'PaperPositionMode','auto')
%print('-dpsc2',strcat('fig_usurf.eps'))

figure(3)
imagesc(x/1000,flipud(-y)/1000,flipud(outline))
title('outline')
xlabel('x (km)'), ylabel('y (km)')
%set(gcf,'PaperPositionMode','auto')
%print('-dpsc2',strcat('fig_outline.eps'))

spera = 31556926.0;
Phi0arr = zeros(size(icemask));
Phi = zeros(size(icemask));
for i=1:length(x)
    for j=1:length(y)
        Phi0arr(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/spera; 
        if outline(i,j)>0.5, Phi(i,j) = Phi0arr(i,j); end
        if usurf(i,j)==0.0, usurf(i,j)=topg(i,j); end
    end
end

figure(4)
imagesc(x/1000,flipud(-y)/1000,flipud(Phi*spera)), colorbar
title('water input  (m/a)')
xlabel('x (km)'), ylabel('y (km)')
%set(gcf,'PaperPositionMode','auto')
%print('-dpsc2',strcat('fig_melt.eps'))

return

W0 = zeros(size(topg));
ts = 0.0;
te = tyears*spera;
W = conservewater(x,y,topg,usurf,outline,W0,Phi,ts,te,5);

