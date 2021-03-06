%% Figure5

% Longwave/Shortwave radiative cooling as a function
% of the column water vapor at fixed surf. temp. Ts
% 1. At varying surface mixing ratio
% 2. At varying shape

%% 0. Constants and parameters

% 0.1 Constants

cp      =   1005;           % Specific heat capacity of dry air at constant pressure (J/kg/K)
D       =   1.66;           % Schwarzschild two-stream diffusivity factor
g       =   9.81;           % Surface gravity constant (m/s^2)
Gamma_d =   g/cp;           % Dry adiabatic lapse rate (K/m)
kLW     =   0.1;            % LW absorbtion coefficient of water vapor (kg^-1.m^2)
kSW     =   0.01;           % SW absorbtion coefficient of water vapor (kg^-1.m^2)
Lv      =   2.5e6;          % Latent heat of vaporization of water vapor (J/kg)
mu      =   pi/4;           % Insolation-weighted zenith angle at the Equator
n_RCE   =   3;              % Mixing ratio profile exponent based on Chuuk Lagoon
ps      =   1000;           % Surface pressure (hPa)
sig     =   5.67e-8;        % Stefan-Boltzmann constant (W/m^2/K^4)
S       =   sig*280^4;      % Insolation for an effective emission temperature of 280K (W/m^2)
Ts      =   300;            % Surface temperature (K)
Gm      =   Gamma_m(Ts,ps); % Average moist adiabatic lapse rate (K/m)
Rd      =   287;            % Dry specific gas constant (J/kg/K)
rs_RCE  =   r_sat(ps,Ts);   % Saturation mixing ratio 

% 0.2 Parameters

alpha0  =4*Gm*Rd/g;             % Exponent relating temperature and optical depth
alpha_RCE=alpha0/(n_RCE+2);     % RCE value of this exponent
rh      =linspace(0.01,60,100); % Column water vapor
Lrh     =length(rh);
eps     =kSW/(kLW*mu*D);        % Ratio of shortwave to longwave optical depth

%% 1. Varying surface mixing ratio, fixed shape

tau_fixn=D*(n_RCE+1)*kLW/(n_RCE+2)*rh; % Optical depth of atm at fixed n
QS_fixn=S*(rh.^0-exp(-eps*tau_fixn)); % SW atm radiative cooling

QL_fixn=zeros(Lrh,1); % Longwave atm cooling

for irh=1:Lrh
    
    ta=tau_fixn(irh);
    y1 = @(x) (x/ta).^alpha_RCE.*exp(x-ta);
    y2 = @(x) (x/ta).^alpha_RCE.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    QL_fixn(irh)=sig*Ts^4*(exp(-ta)-1+Y1+Y2);
    
end

%% 2. Varying shape, fixed surface mixing ratio

rsh=rs_RCE*100*ps/g; % Column water vapor if r(p)=rs
tau_fixrs=kLW*D*rsh*rh./(rsh+rh); % Optical depth of atm at fixed rs
alp_fixrs=alpha0*rh./(rsh+rh); % Exponent alpha at fixed rs
QS_fixrs=S*(rh.^0-exp(-eps*tau_fixrs)); % SW atm radiative cooling

QL_fixrs=zeros(Lrh,1); % LW atm radiative cooling

for irh=1:Lrh
    
    ta=tau_fixrs(irh);
    al=alp_fixrs(irh);
    y1 = @(x) (x/ta).^al.*exp(x-ta);
    y2 = @(x) (x/ta).^al.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    QL_fixrs(irh)=sig*Ts^4*(Y1+Y2+exp(-ta)-1);
    
end

Q_fixn = QL_fixn-QS_fixn';
Q_fixrs = QL_fixrs-QS_fixrs';

%% 3. Plot

figure('position',[100 100 700 300]);

subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
plot(rh,Q_fixn,'Linewidth',1,'color',[0 0 0]);
plot(rh,QL_fixn,'Linewidth',1,'Linestyle','--','color',[0 0 0]);
[~,ind_max]=max(QL_fixn);
plot(rh(ind_max),QL_fixn(ind_max),'+','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixn);
plot(rh(ind_max),Q_fixn(ind_max),'+','Linewidth',1,'color',[0 0 0]);
box on;
ylim([0 500]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
ylabel('Radiative cooling, $\hat{Q}$ [W m$^{-2}$]','Interpreter','Latex');
title('a) Varied $r_s$','Fontweight','normal','Interpreter','Latex');
LEG=legend('$\hat{Q}=\hat{Q}_L+\hat{Q}_S$','$\hat{Q}_L$');
set(LEG,'Fontsize',11,'Location','SouthEast','Interpreter','Latex');

subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
plot(rh,Q_fixrs,'Linewidth',1,'color',[0 0 0]);
plot(rh,QL_fixrs,'Linewidth',1,'Linestyle','--','color',[0 0 0]);
[~,ind_max]=max(QL_fixrs);
plot(rh(ind_max),QL_fixrs(ind_max),'+','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixrs);
plot(rh(ind_max),Q_fixrs(ind_max),'+','Linewidth',1,'color',[0 0 0]);
box on;
ylim([0 500]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
title('b) Varied $n$','Fontweight','normal','Interpreter','Latex');
%LEG=legend('Q_L+Q_S','Q_L');
%set(LEG,'Fontsize',16,'Location','SouthEast');

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure05_LW_SW_gray.pdf']);