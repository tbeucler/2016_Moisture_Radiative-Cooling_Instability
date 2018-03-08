%% Figure 14

% Longwave radiative cooling, surface longwave flux and outgoing longwave
% radiation, as a funct. of the column water vapor at fixed surf. temp. Ts
% 1. At varying surface mixing ratio
% 2. At varying shape

%% 0. Constants and parameters

% 0.1 Constants

cp      =   1005;           % Specific heat capacity of dry air at constant pressure (J/kg/K)
D       =   1.66;           % Schwarzschild two-stream diffusivity factor
g       =   9.81;           % Surface gravity constant (m/s^2)
Gamma_d =   g/cp;           % Dry adiabatic lapse rate (K/m)
k1      =   1.0;            % Absorbtion coefficient of water vapor in bands (kg^-1.m^2)
k2      =   0.02;           % Absorbtion coefficient of water vapor in windows (kg^-1.m^2)
b1      =   0.5;            % Fraction of spectrum in bands
b2      =   1-b1;           % Fraction of spectrum in windows
Lv      =   2.5e6;          % Latent heat of vaporization of water vapor (J/kg)
n_RCE   =   3;              % Mixing ratio profile exponent based on Chuuk Lagoon
ps      =   1000;           % Surface pressure (hPa)
sig     =   5.67e-8;        % Stefan-Boltzmann constant (W/m^2/K^4)
Ts      =   300;            % Surface temperature (K)
Gm      =   Gamma_m(Ts,ps); % Average moist adiabatic lapse rate
Rd      =   287;            % Dry specific gas constant (J/kg/K)
rs_RCE  =   r_sat(ps,Ts);   % Saturation mixing ratio

% 0.2 Parameters

alpha0  =   4*Gm*Rd/g;          % Exponent relating temperature and optical depth
alpha_RCE=  alpha0/(n_RCE+2);   % RCE value of this exponent
rh      =linspace(0.01,95,100); % Column water vapor
Lrh     =length(rh);

%% 1. Varying surface mixing ratio, fixed shape

tau_fixn1=D*(n_RCE+1)*k1/(n_RCE+2)*rh; % Optical depth of atm at fixed n in bands
tau_fixn2=D*(n_RCE+1)*k2/(n_RCE+2)*rh; % Optical depth of atm at fixed n in windows

SLW_fixn=zeros(Lrh,1); % Surface net longwave flux at fixed n
OLR_fixn=zeros(Lrh,1); % Outgoing longwave radiation at fixed n

for irh=1:Lrh
    
    % bands
    ta=tau_fixn1(irh);
    y1 = @(x) (x/ta).^alpha_RCE.*exp(x-ta);
    y2 = @(x) (x/ta).^alpha_RCE.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    SLW_fixn(irh)=b1*sig*Ts^4*(1-Y1);
    OLR_fixn(irh)=b1*sig*Ts^4*(exp(-ta)+Y2);
    
    % windows
    ta=tau_fixn2(irh);
    y1 = @(x) (x/ta).^alpha_RCE.*exp(x-ta);
    y2 = @(x) (x/ta).^alpha_RCE.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    SLW_fixn(irh)=SLW_fixn(irh)+b2*sig*Ts^4*(1-Y1);
    OLR_fixn(irh)=OLR_fixn(irh)+b2*sig*Ts^4*(exp(-ta)+Y2);
    
end

Q_fixn=OLR_fixn-SLW_fixn; % Radiative cooling at fixed n

%% 2. Varying shape, fixed surface mixing ratio

rsh=rs_RCE*100*ps/g; % Column water vapor if r(p)=rs
tau_fixrs1=k1*D*rsh*rh./(rsh+rh); % Optical depth of atm at fixed rs in bands
alp_fixrs1=alpha0*rh./(rsh+rh); % Exponent alpha at fixed rs in bands

tau_fixrs2=k2*D*rsh*rh./(rsh+rh); % Optical depth of atm at fixed rs in windows
alp_fixrs2=alpha0*rh./(rsh+rh); % Exponent alpha at fixed rs in windows

SLW_fixrs=zeros(Lrh,1); % Surface net longwave flux at fixed rs
OLR_fixrs=zeros(Lrh,1); % Outgoing longwave radiation at fixed rs

for irh=1:Lrh
    
    % bands
    ta=tau_fixrs1(irh);
    al=alp_fixrs1(irh);
    y1 = @(x) (x/ta).^al.*exp(x-ta);
    y2 = @(x) (x/ta).^al.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    SLW_fixrs(irh)=b1*sig*Ts^4*(1-Y1);
    OLR_fixrs(irh)=b1*sig*Ts^4*(exp(-ta)+Y2);
    
    % windows
    ta=tau_fixrs2(irh);
    al=alp_fixrs2(irh);
    y1 = @(x) (x/ta).^al.*exp(x-ta);
    y2 = @(x) (x/ta).^al.*exp(-x);
    Y1=integral(y1,0,ta);
    Y2=integral(y2,0,ta);
    SLW_fixrs(irh)=SLW_fixrs(irh)+b2*sig*Ts^4*(1-Y1);
    OLR_fixrs(irh)=OLR_fixrs(irh)+b2*sig*Ts^4*(exp(-ta)+Y2);
    
end

Q_fixrs=OLR_fixrs-SLW_fixrs; % Radiative cooling at fixed rs

%% 3. Plot

figure('position',[100 100 700 300]);

subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
plot(rh,SLW_fixn,'Linewidth',1,'color',[0 0 1]);
plot(rh,OLR_fixn,'Linewidth',1,'color',[1 0 0]);
plot(rh,Q_fixn,'--','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixn);
plot(rh(ind_max),Q_fixn(ind_max),'+','Linewidth',1,'color',[0 0 0]);
box on;
ylim([0 500]);
xlim([0 90]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]', 'Interpreter', 'Latex');
ylabel('Radiative flux [W m$^{-2}$]', 'Interpreter', 'Latex');
title('a) Varied $r_s$','FontWeight','normal', 'Interpreter', 'Latex');
LEG=legend('SLW','OLR','$\hat{Q}_L$');
set(LEG,'Position',[0.31    0.65    0.14    0.23], 'Fontsize',11, 'Interpreter', 'Latex');

subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
plot(rh,SLW_fixrs,'Linewidth',1,'color',[0 0 1]);
plot(rh,OLR_fixrs,'Linewidth',1,'color',[1 0 0]);
plot(rh,Q_fixrs,'--','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixrs);
plot(rh(ind_max),Q_fixrs(ind_max),'+','Linewidth',1,'color',[0 0 0]);
box on;
ylim([0 500]);
xlim([0 90]);
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]', 'Interpreter', 'Latex');
title('b) Varied $n$','FontWeight','normal', 'Interpreter', 'Latex');
%LEG=legend('SLW','OLR','Q_L');
%set(LEG,'Location','SouthEast','Fontsize',10);

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure14_LW_2bandgray.pdf']);