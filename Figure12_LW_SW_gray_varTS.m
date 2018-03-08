%% Figure12

% Longwave/Shortwave radiative cooling as a function
% of the column water vapor for different surface temperatures
% 1. At varying surface mixing ratio
% 2. At varying shape

%% 0. Constants and parameters

% 0.1 Constants

cp      =   1005;               % Specific heat capacity of dry air at constant pressure (J/kg/K)
D       =   1.66;               % Schwarzschild two-stream diffusivity factor 
g       =   9.81;               % Surface gravity constant (m/s^2) 
Gamma_d =   g/cp;               % Dry adiabatic lapse rate (K/m)
kLW     =   0.1;                % LW absorbtion coefficient of water vapor (kg^-1.m^2)
kSW     =   0.01;               % SW absorbtion coefficient of water vapor (kg^-1.m^2)
Lv      =   2.5e6;              % Latent heat of vaporization of water vapor (J/kg)
mu      =   pi/4;               % Insolation-weighted zenith angle at the Equator
n_RCE   =   3;                  % Mixing ratio profile exponent based on Chuuk Lagoon
ps      =   [1000 1000 1000];   % Surface pressure (hPa)
sig     =   5.67e-8;            % Stefan-Boltzmann constant (W/m^2/K^4)
S       =   sig*280^4;          % Insolation for an effective emission temperature of 280K (W/m^2)
Ts      =   [285 300 315];      % Surface temperature (K)
LTs     =   length(Ts);         % Number of surface temperatures
Gm      =   Gamma_m(Ts,ps);     % Average moist adiabatic lapse rate (K/m)
Rd      =   287;                % Dry specific gas constant (J/kg/K)
rs_RCE  =   r_sat(ps,Ts);       % Saturation mixing ratio

% 0.2 Parameters

alpha0      =   4*Gm*Rd/g;              % Exponent relating temperature and optical depth
alpha_RCE   =   alpha0/(n_RCE+2);       % RCE value of this exponent
rh          =   linspace(0.01,60,100);  % Column water vapor (kg/m^2)
Lrh         =   length(rh);
eps         =   kSW/(kLW*mu*D);         % Ratio of shortwave to longwave optical depth

%% 1. Varying surface mixing ratio, fixed shape

tau_fixn=D*(n_RCE+1)*kLW/(n_RCE+2)*rh; % Optical depth of atm at fixed n

QS_fixn=zeros(Lrh,LTs); % Shortwave atm cooling
QL_fixn=zeros(Lrh,LTs); % Longwave atm cooling

for irh=1:Lrh
    
    ta=tau_fixn(irh);
    
    for iTs=1:LTs
        
        alp=alpha_RCE(iTs);
        
        y1 = @(x) (x/ta).^alp.*exp(x-ta);
        y2 = @(x) (x/ta).^alp.*exp(-x);
        Y1=integral(y1,0,ta);
        Y2=integral(y2,0,ta);
        QL_fixn(irh,iTs)=sig*Ts(iTs)^4*(exp(-ta)-1+Y1+Y2);
        QS_fixn(irh,iTs)=S*(1-exp(-eps*ta)); % SW atm radiative cooling
        
    end
    
end

%% 2. Varying shape, fixed surface mixing ratio

%rsh=rs_RCE(2)*100*ps(2)/g; % Column water vapor if r(p)=rs
%tau_fixrs=kLW*D*rsh*rh./(rsh+rh); % Optical depth of atm at fixed rs

QS_fixrs=zeros(Lrh,LTs); % SW atm radiative cooling
QL_fixrs=zeros(Lrh,LTs); % LW atm radiative cooling

for irh=1:Lrh
    
    for iTs=1:LTs
        
        rsh=rs_RCE(iTs)*100*ps(iTs)/g; % Column water vapor if r(p)=rs
        ta=kLW*D*rsh*rh(irh)./(rsh+rh(irh)); % Optical depth of atm at fixed rs
% twc - commented these lines out and moved some to before loop -- I don't
% think we want to be using varied surface mixing ratio at different
% temperatures, do we??
%
%        ta = tau_fixrs(irh);
        alp=alpha0(iTs)*rh(irh)./(rsh+rh(irh)); % Exponent alpha at fixed rs
        QS_fixrs(irh,iTs)=S*(1-exp(-eps*ta)); % SW atm radiative cooling
        
        
        y1 = @(x) (x/ta).^alp.*exp(x-ta);
        y2 = @(x) (x/ta).^alp.*exp(-x);
        Y1=integral(y1,0,ta);
        Y2=integral(y2,0,ta);
        QL_fixrs(irh,iTs)=sig*Ts(iTs)^4*(Y1+Y2+exp(-ta)-1);
        
    end
    
end

Q_fixn = QL_fixn-QS_fixn;
Q_fixrs = QL_fixrs-QS_fixrs;

%% 3. Plot

figure('position',[100 100 700 300]);
cmap                = [0.3 0.3 0.6; 0 0 0; 0.6 0.3 0.3];
tau_labels          = {'1';'2';'3';'4';'5';'6';'7'};
rh_taulabel_fixn    = interp1(tau_fixn,rh,(1:7));
rh_taulabel_fixrs   = interp1(tau_fixrs,rh,(1:7));

subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
HP1 = plot(rh,Q_fixn,'Linewidth',1);
set(HP1, {'color'}, num2cell(cmap, 2));
[~,ind_max]=max(Q_fixn);
for iTs=1:LTs
    plot(rh(ind_max(iTs)),Q_fixn(ind_max(iTs),iTs),'+','Linewidth',1,'color',cmap(iTs,:));
end
box on;
% text(0,485,'\tau_s=','Color',[0 0.5 0],'HorizontalAlignment','left','VerticalAlignment','top');
% for i=1:length(tau_labels)
%     text(rh_taulabel_fixn(i),485,char(tau_labels(i)),'Color',[0 0.5 0],'HorizontalAlignment','center','VerticalAlignment','top');
%     line([rh_taulabel_fixn(i) rh_taulabel_fixn(i)],[0 10],'Color',[0 0.5 0],'LineWidth',1);
%     line([rh_taulabel_fixn(i) rh_taulabel_fixn(i)],[490 500],'Color',[0 0.5 0],'LineWidth',1);
% end
ylim([0 500]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','latex');
ylabel('Radiative cooling, $\hat{Q}$ [W m$^{-2}$]','Interpreter','latex');
title('a) Varied $r_s$','Interpreter','latex','Fontweight','normal');
LEG=legend([num2str(Ts(1)),'K'],[num2str(Ts(2)),'K'],[num2str(Ts(3)),'K']);
set(LEG,'Fontsize',11,'Position',[0.175  0.17  0.145  0.185]);

subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
HP2 = plot(rh,Q_fixrs,'Linewidth',1);
set(HP2, {'color'}, num2cell(cmap, 2));
[~,ind_max]=max(Q_fixrs);
for iTs=1:LTs
    plot(rh(ind_max(iTs)),Q_fixrs(ind_max(iTs),iTs),'+','Linewidth',1,'color',cmap(iTs,:));
end
box on;
% text(0,485,'\tau_s=','Color',[0 0.5 0],'HorizontalAlignment','left','VerticalAlignment','top');
% for i=1:length(tau_labels)
%     text(rh_taulabel_fixrs(i),485,char(tau_labels(i)),'Color',[0 0.5 0],'HorizontalAlignment','center','VerticalAlignment','top');
%     line([rh_taulabel_fixrs(i) rh_taulabel_fixrs(i)],[0 10],'Color',[0 0.5 0],'LineWidth',1);
%     line([rh_taulabel_fixrs(i) rh_taulabel_fixrs(i)],[490 500],'Color',[0 0.5 0],'LineWidth',1);
% end
ylim([0 500]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','latex');
title('b) Varied $n$','Interpreter','latex','Fontweight','normal');
%LEG=legend([num2str(Ts(1)),'K'],[num2str(Ts(2)),'K'],[num2str(Ts(3)),'K']);
%set(LEG,'Fontsize',16,'Location','SouthEast');

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure12_LW_SW_gray_varTS.pdf']);