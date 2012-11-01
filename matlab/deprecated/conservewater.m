function W = conservewater(x,y,b,h,outline,W0,Phi,ts,te,Nmin)
% CONSERVE Runs the minimal advection-diffusion model for subglacial
% hydrology from time  ts  to time  te  using at least  Nmin  time steps:
%   W_t + div q = Phi
%   q = - (K / rhow g) (grad psi) W
%   psi = P + rhow g (b + W)
%   P = s P_i = s rhoi g H
% where 
%   b = bed elevation
%   h = surface elevation
%   H = h - b  if  floatmask==0
% Note that
%   v = (alpha,beta) = - s K r grad h - K (1 - s r) grad b
% and the flux  q  is a sum of advective and diffusive components:
%   q = v W - K W grad W
% so the equation is
%   W_t + div (v W) = div (K W grad W) + Phi
% If floatmask==1 then W=0 at that point.  (FIXME: account for water lost.)
% Form:  W = conserve(x,y,b,h,floatmask,W0,Phi,ts,te,Nmin)

%% OUTPUT FOLDER
foldername = 'OUT_test_1';
mkdir(foldername);


%% PARAMETERS
% constants
rhoi = 910.0;       % ice density (kg m-3)
rhow = 1028.0;      % water density (kg m-3)
spera = 31556926.0; % seconds per yr
r = rhoi / rhow;
g = 9.81;           % gravitation acceleration (m s-2)

% model parameters
Kmin = 5e-3;        % low conductivity (m s-1)
Kmax = 5e-2;        % high conductivity (m s-1)
A = 3.1689d-24;     % viscosity parameter (Pa-3 s-1)
c1 = 0.500;         % cavitation constant (m-1)
c2 = 0.040;         % creep closure constant (dimensionless)
vb = 50.0/spera;    % ice velocity (m a-1)
Ymin = 0.0001;      % minimum W required in Pw formulation (avoid division by zero)
Wr = 1.0;           % basal roughness scale (W>Wr: sheet instability --> channelization)

%% ARRAY DECLARATION
Qadv = zeros(size(b));      % advective flux
Qdiff = zeros(size(b));     % diffusive flux
scale = zeros(size(b));     % water pressure / overburden * r
bwp = zeros(size(b));       % basal water pressure
MR = zeros(size(b));        % melt rate
alpha = zeros(size(b));     % factor in flux divergence (x) 
beta = zeros(size(b));      % factor in flux divergence (y)
K = zeros(size(b));         % conductivity
Wt = zeros(size(b));        % rate of change of water thickness
pressure = zeros(size(b));  % overpressure (2), normal pressure (1), underpressure (0)
E = zeros(size(b));         % non-filled space of cavity
Et = zeros(size(b));		% E/dt


[Mx, My] = size(W0);
dx = x(2) - x(1);
dy = y(2) - y(1);

%% BED & SURFACE SLOPE
dhdx  = (h(2:end,:)-h(1:end-1,:)) / dx;
dbdx  = (b(2:end,:)-b(1:end-1,:)) / dx;
dhdy  = (h(:,2:end)-h(:,1:end-1)) / dy;
dbdy  = (b(:,2:end)-b(:,1:end-1)) / dy;

fprintf('running ...\n')
fprintf('                                ')
fprintf('max W (m)  av W (m)  vol(10^6 km^3)  rel-bal\n')

t = ts;
W = W0;
volW = 0.0;
dA = dx*dy;

count = 1;
counter = 1;

while t<te
  dtmax = (te - ts) / Nmin;
  
  %% DETERMINE WATER PRESSURE
  for i=2:length(x)-1
    for j=2:length(y)-1
      Wij = W(i,j);
      Wtij = Wt(i,j);
      Eij = E(i,j);
      Etij = Et(i,j);
      
      Yij = Eij + Wij;
      Ytij = Wtij - Etij;
      
      Po = rhoi*g*(h(i,j)-b(i,j));
      if (h(i,j)==0.0), Po = 0.0; end 
      Vcav = c1*abs(vb)*(Wr-Yij);
      if Vcav<0.0, Vcav = 0.0; end     
      
      % OVERPRESSURE: 
      % - occurs when dYdt>Vcav or W>Wr
      % - K = Kmax when W>Wr, K = Kmin when W<=Wr
      % - Pw = Po
      % UNDERPRESSURE:
      % - occurs when closure of the system cannot keep up with negative dWdt
      % - K = Kmin
      % - Pw = 0.0
      % NORMAL PRESSURE:
      % - occurs in any other case
      % - K = Kmin
      % - Pw = Po - ((Vcav-Wt)/(c2*A*W))^(1/3), follows from Wt = Vcav - Vcreep
      
      if (Wij>Wr || Ytij>Vcav)
          Pw = Po;
          pressure(i,j) = 2.0;
          if Wij>Wr
              K(i,j) = Kmax; 
          else
              K(i,j) = Kmin;
          end
      elseif ((Vcav-Ytij)/(c2*A*(Yij+Ymin)))^(1/3)>Po;
          Pw = 0.0;
          pressure(i,j) = 0.0;
          K(i,j) = Kmin;
      else
          Pw = Po - ((Vcav-Ytij)/(c2*A*(Yij+Ymin)))^(1/3);
          pressure(i,j) = 1.0;
          K(i,j) = Kmin;
      end
      
      if Po>0.0
        s = Pw/Po;
      else
        s = 0.0;
      end
      
      % BOUNDARY CONDITION (Pw = 0.0)
      if (outline(i,j)<0.5)
          Pw = 0.0;
          s = 0.0;
      end
      
      scale(i,j) = s*r;
      bwp(i,j) = Pw;
    end
  end
  
  
  %% EVOLVE WATER THICKNESS
  
  % staggered s, K and W
  Sea = 0.5 * (scale(2:end,:) + scale(1:end-1,:));
  Sno = 0.5 * (scale(:,2:end) + scale(:,1:end-1));
  Kea = 0.5 * (K(2:end,:) + K(1:end-1,:));
  Kno = 0.5 * (K(:,2:end) + K(:,1:end-1));
  Wea = 0.5 * (W(2:end,:) + W(1:end-1,:));
  Wno = 0.5 * (W(:,2:end) + W(:,1:end-1));
  
  % alpha, beta and melt rate (MR)
  for i=2:length(x)-1
      for j=2:length(y)-1
        alpha(i,j) = - Kea(i,j) * (Sea(i,j) * dhdx(i,j) + (1.0-Sea(i,j)) * dbdx(i,j));
        beta(i,j)  = - Kno(i,j) * (Sno(i,j) * dhdy(i,j) + (1.0-Sno(i,j)) * dbdy(i,j));

        % sinusoidal summer melt peak
        % MR(i,j) = -2.0*pi*Phi(i,j)*sin(4.0*pi*t/spera-pi/2.0);
        if (MR(i,j)<0.0 || mod(t/spera,1.0)<0.25 || mod(t/spera,1.0)>0.75), MR(i,j)=0.0; end

        % discrete melt rate
        if (mod(t/spera,1.0)>0.25 && mod(t/spera,1.0)<0.5), MR(i,j) = 3.0*Phi(i,j); end  
        if (mod(t/spera,1.0)>0.5 && mod(t/spera,1.0)<0.75), MR(i,j) = 1.0*Phi(i,j); end
      end
  end
  
  % determine time-step (satisfying criteria for advection and diffusion)
  dtCFL = 0.5 / (max(max(alpha))/dx + max(max(beta))/dy);
  
  dt = min(te-t,dtmax);

  maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  dtDIFF = 0.25 / (Kmax * maxW * (1.0/dx^2 + 1.0/dy^2));
  
  dt = min([dt dtCFL dtDIFF]);

  fprintf('  t = %8.3f a [dt = %6.3f]:  ',(t+dt)/spera,dt/spera)

  
  
  nux = dt / dx;        nuy = dt / dy;
  Wnew = W;  % copies unaltered zero b.c.s
  
  inputvol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      Wij = W(i,j);
      Eij = E(i,j);
      Wtij = Wt(i,j);
      Yij = Wij + Eij;
      
      % update empty space (nonzero in case of underpressure)
      Po = rhoi*g*(h(i,j)-b(i,j));
      if (pressure(i,j)==0.0 && outline(i,j)>0.5)
          E(i,j) = Eij + (c1*(Wr-Yij)*abs(vb) - c2*A*Po^3*Yij - Wtij) * dt;
      else
          E(i,j) = 0.0;
      end
      Et(i,j) = E(i,j)/dt;
      
      
      % update water thickness following mass conservation
      mux = K(i,j) * dt / dx^2;  muy = K(i,j) * dt / dy^2;  
   
      upe = up(alpha(i,j),  Wij,     W(i+1,j));
      upw = up(alpha(i-1,j),W(i-1,j),Wij);
      upn = up(beta(i,j),   Wij,     W(i,j+1));
      ups = up(beta(i,j-1), W(i,j-1),Wij);
      inputdepth = dt * MR(i,j);
      
      Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + ...
          mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
          muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1))) + ...
          inputdepth;
      
      % rate of change in water thickness
      Wt(i,j) = (Wnew(i,j)-Wij)/dt;
      
      % advective flux
      Qadvx = abs(upe+upw)/2.0;
      Qadvy = abs(upn+ups)/2.0;
      Qadv(i,j) = sqrt(Qadvx^2+Qadvy^2);
      
      % diffusive flux
      Qdiffx = abs(K(i,j)*Wij*((Wea(i,j)-Wea(i-1,j))/dx));
      Qdiffy = abs(K(i,j)*Wij*((Wno(i,j)-Wno(i-1,j))/dy));
      Qdiff(i,j) = sqrt(Qdiffx^2+Qdiffy^2);
      
      inputvol = inputvol + inputdepth * dA;
    end
  end

  losevol = 0.0;
  for i=2:length(x)-1
    for j=2:length(y)-1
      if outline(i,j) < 0.5
        losevol = losevol + Wnew(i,j)*dA; 
        Wnew(i,j) = 0.0;
      end
    end
  end

  sumnew = sum(sum(Wnew));
  fprintf('%7.3f  %7.3f',max(max(Wnew)),sumnew/(Mx*My))

  volnew = sumnew * dA;
  balance = (volW + inputvol - losevol) - volnew;
  fprintf('        %.3f        %.0e\n',...
            volnew/(1e9*1e6),abs(balance/volnew))

  volW = volnew;
  W = Wnew;
  t = t + dt;
  
  Nice = 0;
  for i=1:length(x)
    for j=1:length(y)
        if outline(i,j)>0.5, Nice = Nice + 1; end
        if outline(i,j)<0.5, pressure(i,j) = 3.0; end
    end
  end
  
  Wmean(counter,1) = t/spera;
  Pwmean(counter,1) = t/spera;
  Mmean(counter,1) = t/spera;
  PW_under(counter,1) = t/spera;
  PW_over(counter,1) = t/spera;
  PW_normal(counter,1) = t/spera;
  Emean(counter,1) = t/spera;
  VolWmean(counter,1) = t/spera;
  Vollosemean(counter,1) = t/spera;
  Volinputmean(counter,1) = t/spera;
  Wmean(counter,2) = sum(sum(W(outline>0.5)))/Nice;
  Pwmean(counter,2) = sum(sum(bwp(outline>0.5)))/Nice;
  Mmean(counter,2) = sum(sum(MR(outline>0.5)))/Nice;
  Emean(counter,2) = sum(sum(E(outline>0.5)))/Nice;
  VolWmean(counter,2) = volW;
  Vollosemean(counter,2) = losevol;
  Volinputmean(counter,2) = inputvol;
  PW_under(counter,2) = length(find(pressure==0.0))/Nice;
  PW_over(counter,2) = length(find(pressure==2.0))/Nice;
  PW_normal(counter,2) = length(find(pressure==1.0))/Nice;
  
  counter = counter + 1;
  
  snapperiod = 0.05;
  if (t/spera>4.0 && mod(t/(snapperiod*spera),1)<0.5 && mod((t-dt)/(snapperiod*spera),1)>0.5)
%       figure
%       imagesc(x/1000,flipud(-y)/1000,flipud(W)), colorbar
%       title(sprintf('Water amount from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_W_',int2str(20*snapperiod*count),'.eps'))
% 
% 
%       imagesc(x/1000,flipud(-y)/1000,flipud(Qadv)), colorbar
%       title(sprintf('Advection from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_Qadv_',int2str(20*snapperiod*count),'.eps'))
% 
% 
%       imagesc(x/1000,flipud(-y)/1000,flipud(Qdiff)), colorbar
%       title(sprintf('Diffusion from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_Qdif_',int2str(20*snapperiod*count),'.eps'))
% 
% 
%       imagesc(x/1000,flipud(-y)/1000,flipud(E)), colorbar
%       title(sprintf('Empty space from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_E_',int2str(20*snapperiod*count),'.eps'))
% 
% 
%       imagesc(x/1000,flipud(-y)/1000,flipud(bwp)), colorbar
%       title(sprintf('Water pressure from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_bwp_',int2str(20*snapperiod*count),'.eps'))
%       
% %       imagesc(x/1000,flipud(-y)/1000,flipud(alpha)), colorbar
% %       title(sprintf('ALPHA from CONSERVE after time of %.3f year',(t)/spera))
% %       xlabel('x (km)'), ylabel('y (km)')
% %       set(gcf,'PaperPositionMode','auto')
% %       print('-dpsc2',strcat(foldername,'/fig_test_alpha_',int2str(20*snapperiod*count),'.eps'))
% %       
% %       imagesc(x/1000,flipud(-y)/1000,flipud(beta)), colorbar
% %       title(sprintf('BETA from CONSERVE after time of %.3f year',(t)/spera))
% %       xlabel('x (km)'), ylabel('y (km)')
% %       set(gcf,'PaperPositionMode','auto')
% %       print('-dpsc2',strcat(foldername,'/fig_test_beta_',int2str(20*snapperiod*count),'.eps'))
%       
%       imagesc(x/1000,flipud(-y)/1000,flipud(pressure)), colorbar
%       title(sprintf('Pressure regime from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_pressure_',int2str(20*snapperiod*count),'.eps'))
%       
%       imagesc(x/1000,flipud(-y)/1000,flipud(K)), colorbar
%       title(sprintf('K from CONSERVE after time of %.3f year',(t)/spera))
%       xlabel('x (km)'), ylabel('y (km)')
%       set(gcf,'PaperPositionMode','auto')
%       print('-dpsc2',strcat(foldername,'/fig_test_K_',int2str(20*snapperiod*count),'.eps'))
%       close(gcf);
      count = count+1;
  end
end


A = dataset({Wmean(:,1),'t'},{Wmean(:,2),'W'},{Pwmean(:,2),'Pw'},{Mmean(:,2),'M'},{PW_under(:,2),'PW_under'},{PW_over(:,2),'PW_over'},{PW_normal(:,2),'PW_normal'},{Emean(:,2),'E'},{VolWmean(:,2),'VolW'},{Vollosemean(:,2),'VolLose'},{Volinputmean(:,2),'VolInput'});
export(A,'file',strcat(foldername,'/out_text.txt'));
  
   function z = up(v,L,R)
   if v >= 0
     z = v * L;
   else
     z = v * R;
   end
