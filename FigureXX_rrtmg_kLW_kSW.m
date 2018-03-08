%% FigureXX_rrtmg_kLW_kSW
% script to calculate LW and SW effective absorption coefficients as a
% function of column water vapor
% based on the slope of -log(transmissivity) versus r-hat

thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

datadir = [basedir(1:strfind(basedir,'Figure_scripts')-1) 'rrtmg_idealized/'];

% load precipitable water arrays for two variants (var-rs and var-n)
load([datadir 'pw_idealized_chuuk.mat']);
rhat = pw(:,1);

n = 3;      % pressure-water vapor scaling exponent
D = 5/3;    % diffusivity for LW radiation
mu = pi/4;  % cosine of zenith angle for SW radiation

% load sw fluxes for var-rs case
fluxes_var_rs = load([datadir 'output_summary_lw-sw_var_rs.txt']);
swnt = fluxes_var_rs(:,4);
swns = fluxes_var_rs(:,5);

% load lw fluxes for var-rs case, for both control and TS+1K cases
fluxes_var_rs = load([datadir 'output_summary_lw-transmissivity_var_rs.txt']);
lwnt0 = fluxes_var_rs(:,2);
lwns0 = fluxes_var_rs(:,3);
lwnt1 = fluxes_var_rs(:,4);
lwns1 = fluxes_var_rs(:,5);

% SW transmissivity: simple ratio of net surface flux to net TOA flux
T_SW = swns./swnt;

% LW transmissivity: more complicated, 
% based on ratio of change in TOA flux to change in surface flux for 1 K
% surface warming (with no changes to atmosphere)
T_LW = (lwnt1-lwnt0)./(lwns1-lwns0);

rhat_avg = 0.5*(rhat(1:end-1) + rhat(2:end));

k_SW = (mu*(n+2)/(n+1))*diff(-log(T_SW))./diff(rhat);
k_LW = (1/D*(n+2)/(n+1))*diff(-log(T_LW))./diff(rhat);

%% make plot

figure('Position',[100 100 500 300]);
hold on;
set(gca,'Fontsize',11);
plot(rhat_avg,k_LW,'LineWidth',1);
plot(rhat_avg,k_SW,'LineWidth',1);
xlabel('Column water vapor, $\hat{r}$, [kg m$^{-2}$]','Interpreter','latex');
ylabel('Absorption coefficient, $k$, [m$^{2}$ kg$^{-1}$]','Interpreter','latex');
title('Whole-band absorption coefficients fit using RRTMG','Interpreter','latex','FontWeight','normal');
LEG = legend('$k_{LW}$','$k_{SW}$');
set(LEG,'interpreter','latex')
box on;
xlim([0 90]);

gcfsavepdf([basedir 'FigureXX_rrtmg_kLW_kSW.pdf']);
 
