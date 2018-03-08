%% Figure 13

% The goal is to prove that a threshold in column water vapor corresponds
% to a threshhold in RCE surface temperature
% The RCE column water vapor is loaded from the one-layer toy model of
% Beucler and Emanuel (2016)
% This quantity is then compared to the column water vapor threshold of the
% MRCI at fixed n and fixed r_s

%% 0. Constants and parameters

% 0.1 Constants

cp      =   1005;       % Specific heat capacity of dry air at constant pressure, J/kg/K
D       =   1.66;       % Schwarzschild two-stream diffusivity factor
g       =   9.81;       % Gravity constant, m/s^2
Gamma_d =   g/cp;       % Dry adiabatic lapse rate, K/m
kLW     =   0.1;        % LW absorbtion coefficient of water vapor (kg^-1.m^2)
kSW     =   0.01;       % SW absorbtion coefficient of water vapor (kg^-1.m^2)
Lv      =   2.5e6;      % Latent heat of vaporization of water vapor, J/kg/K
ps      =   1000;       % surface pressure, hPa
mu      =   pi/4;       % Insolation-weighted zenith angle at the Equator
n_RCE   =   3;          % Mixing ratio profile exponent based on Chuuk Lagoon
sig     =   5.67e-8;    % Stefan-Boltzmann constant, W/m^2/K^4
S       =   sig*280^4;  % Insolation for an effective emission temperature of 280K
Rd      =   287;        % Dry specific gas constant, J/kg/K

% 0.2 RCE column water vapor from Beucler and Emanuel (2016)

load('Figure13_data.mat'); % Loads rRCE [kg.m^-^2] and Ts [C]
LTs     =   length(Ts); % Number of surface temperatures
Ts      =   Ts+273.15; % Surface temperature [K]

rs_RCE  =   zeros(LTs,1); % Saturation mixing ratio
Gm      =   zeros(LTs,1); % Moist adiabatic lapse rate

for iTs=1:LTs
    rs_RCE(iTs) =   r_sat(1000,Ts(iTs)); % Saturation mixing ratio
    Gm(iTs)     =   Gamma_m(Ts(iTs),1000); % Average moist adiabatic lapse rate
end

% 0.3 Parameters

alpha0  =   4*Gm*Rd/g;              % Exponent relating temperature and optical depth
alpha_RCE=  alpha0/(n_RCE+2);       % RCE value of this exponent
rh      =   linspace(0.01,50,100);  % Column water vapor
Lrh     =   length(rh);
eps     =   kSW/(kLW*mu*D);         % Ratio of shortwave to longwave optical depth


%% 1. Varying surface mixing ratio, fixed shape

tau_fixn=D*(n_RCE+1)*kLW/(n_RCE+2)*rh; % Optical depth of atm at fixed n

Q_fixn=zeros(Lrh,LTs); % Total clear-sky atmospheric cooling

for irh=1:Lrh
    
    ta=tau_fixn(irh);
    
    for iTs=1:LTs
        
        alp=alpha_RCE(iTs);
        
        y1 = @(x) (x/ta).^alp.*exp(x-ta);
        y2 = @(x) (x/ta).^alp.*exp(-x);
        Y1=integral(y1,0,ta);
        Y2=integral(y2,0,ta);
        QL=sig*Ts(iTs)^4*(exp(-ta)-1+Y1+Y2); % LW atm rad cooling
        QS=S*(1-exp(-eps*ta)); % SW atm radiative cooling
        Q_fixn(irh,iTs)=QL-QS; % Total clear-sky atm rad cooling
        
    end
    
end

rcrit_fixn=zeros(LTs,1); % Threshold for MRCI at fixed n
for iTs=1:LTs
    [~,imax]=max(squeeze(Q_fixn(:,iTs))); % Threshold = Maximum of rad cooling
    rcrit_fixn(iTs)=rh(imax); % Corresponding column water vapor
end

%% 2. Varying shape, fixed surface mixing ratio

Q_fixrs=zeros(Lrh,LTs); % Total clear-sky atm radiative cooling

for irh=1:Lrh
    
    for iTs=1:LTs
        
        rsh=rs_RCE(iTs)*100*ps/g; % Column water vapor if r(p)=rs
        ta=kLW*D*rsh*rh(irh)./(rsh+rh(irh)); % Optical depth of atm at fixed rs
        alp=alpha0(iTs)*rh(irh)./(rsh+rh(irh)); % Exponent alpha at fixed rs
        QS=S*(1-exp(-eps*ta)); % SW atm radiative cooling
        
        
        y1 = @(x) (x/ta).^alp.*exp(x-ta);
        y2 = @(x) (x/ta).^alp.*exp(-x);
        Y1=integral(y1,0,ta);
        Y2=integral(y2,0,ta);
        QL=sig*Ts(iTs)^4*(Y1+Y2+exp(-ta)-1); % LW atm radiative cooling
        
        Q_fixrs(irh,iTs)=QL-QS; % Total clear-sky atm radiative cooling
        
    end
    
end

rcrit_fixrs=zeros(LTs,1); % Threshold for MRCI at fixed rs
for iTs=1:LTs
    [~,imax]=max(squeeze(Q_fixrs(:,iTs))); % Threshold = Maximum of rad cooling
    rcrit_fixrs(iTs)=rh(imax); % Corresponding column water vapor
end

%% 3. Plot

Ts=Ts-273.15; % Surface temperature [C]

figure('Position',[100 100 400 300])
hold on;
set(gca,'Fontsize',11);
plot(Ts,rRCE,'color',[0 0 0],'Linewidth',1);
plot(Ts,rcrit_fixn,'color',[0 0 1],'Linewidth',1);
plot(Ts,rcrit_fixrs,'color',[1 0 0],'Linewidth',1);
xlabel('RCE surface temperature, $T_s$ [$^{\circ}$C]','Interpreter','latex');
ylabel('Colummn water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','latex');
box on;

LEG=legend('RCE $\hat{r}$','Critical $\hat{r}$ for varied $r_s$','Critical $\hat{r}$ for varied $n$');
set(LEG,'Fontsize',11,'Location','NorthWest','Interpreter','latex');

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure13_T_thresholds.pdf']);