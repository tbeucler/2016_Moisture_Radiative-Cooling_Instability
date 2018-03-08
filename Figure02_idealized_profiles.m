%% Figure 2

% The temperature profile is fitted using an average moist adiabatic lapse
% rate
% The mixing ratio profile is fitted using a power law
% Requires observational temperature and mixing ratio profiles (ideally
% averaged over many years): p[hPa],T[C],r[g/kg]

% Optional temperature and mixing ratio profiles in Chuuk Lagoon
% Ensemble average of the 21 past months of January (1996-2016)

close all;
clearvars;
fclose('all');

%% 0. Constants and observational data

% 0.1 Constants

g   =   9.81;   % Gravity constant
Rd  =   287;    % Specific gas constant of dry air

% 0.2 Observational data

load('Figure2_profiles.mat');
Ts  =   T(1)+273.15;    % Surface temperature [K]
ps  =   p(1);           % Surface pressure
rs  =   r(1);           % Surface mixing ratio

%% 1. Fit temperature profile

Gm      =   Gamma_m(Ts,ps);
fit_T   =   (Ts*(p/ps).^(Gm*Rd/g))-273.15; % Temperature fit [C]

%% 2. Fit mixing ratio profile

n       =   3;              % Adjust this value by eye or using proper algorithm 
% A r2 fit in log-log space gives too much importance to low pressure data
% where there is no water vapor and the uncertainty is greatest

fit_r   =   rs*(p/ps).^n;   % Mixing ratio fit [g/kg]

%% 3. Plot

figure('position',[100 100 900 300]); 

subplot(1,3,1);
hold on;
set(gca,'Fontsize',11);
plot(T,p,'Linewidth',1,'color',[0 0 0]);
plot(fit_T,p,'Linewidth',1,'Linestyle','--','color',[0 0.25 1]);
box on;
set(gca,'Ydir','reverse');
xlabel('Temperature, $T(p)$ [$^{\circ}$C]','Interpreter','Latex');
ylabel('Pressure, $p$ [hPa]','Interpreter','Latex');
title('a) Temperature profile','FontWeight','normal','Interpreter','Latex');
xlim([-100 30]);
L=legend('Observed','Fit');
set(L,'Location','SouthWest','Fontsize',11,'Interpreter','Latex');

subplot(1,3,2);
hold on;
set(gca,'Fontsize',11);
plot(r,p,'Linewidth',1,'color',[0 0 0]);
plot(fit_r,p,'Linewidth',1,'Linestyle','--','color',[0 0.25 1]);
box on;
set(gca,'Ydir','reverse');
xlabel('Mixing ratio, $r(p)$ [g/kg]','Interpreter','Latex');
title('b) Mixing ratio profile (linear)','FontWeight','normal','Interpreter','Latex');
%L=legend('Obs','Fit');
%set(L,'Location','NorthEast','Fontsize',10);

subplot(1,3,3);
hold on;
set(gca,'Fontsize',11);
plot(r,p,'Linewidth',1,'color',[0 0 0]);
plot(fit_r,p,'Linewidth',1,'Linestyle','--','color',[0 0.25 1]);
box on;
set(gca,'Ydir','reverse');
set(gca,'Xscale','log');
xlim([1e-3 100])
xlabel('Mixing ratio, $r(p)$ [g/kg]','Interpreter','Latex');
set(gca,'Xtick',[0.001 0.01 0.1 1 10 100]);
title('c) Mixing ratio profile (log)','FontWeight','normal','Interpreter','Latex');
%L=legend('Obs','Fit');
%set(L,'Location','NorthEast','Fontsize',10);

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
%gcfsavepdf([basedir 'Figure02_idealized_profiles.pdf']);