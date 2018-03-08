%% Fig08_rrtmg_idealized
% script to plot tropospheric radiative cooling against column water vapor
% for rrtm calculations, with extremely idealized water vapor profiles
%     r(p) = r_s (p/p_s)^n
% column water vapor is controlled by varying either r_s or n

thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

datadir = [basedir(1:strfind(basedir,'Figure_scripts')-1) 'rrtmg_idealized/'];

% load precipitable water arrays for two variants (var-rs and var-n)
load([datadir 'pw_idealized_chuuk.mat']);

% load fluxes for var-rs and var-n cases
fluxes_var_rs = load([datadir 'output_summary_lw-sw_var_rs.txt']);
fluxes_var_n = load([datadir 'output_summary_lw-sw_var_n.txt']);

% scale factor for sw fluxes (should depend on zenith angle used and latitude)
sw_scale_factor = 4/pi^2;
% 4/pi^2 appropriate for equatorial equinox with cos(zeta)=pi/4
sw_heating_var_rs = sw_scale_factor*(fluxes_var_rs(:,4)-fluxes_var_rs(:,5));
sw_heating_var_n = sw_scale_factor*(fluxes_var_n(:,4)-fluxes_var_n(:,5));

QL_fixn = fluxes_var_rs(:,2)-fluxes_var_rs(:,3);
QL_fixrs = fluxes_var_n(:,2)-fluxes_var_n(:,3);

Q_fixn = fluxes_var_rs(:,2)-fluxes_var_rs(:,3)-sw_heating_var_rs;
Q_fixrs = fluxes_var_n(:,2)-fluxes_var_n(:,3)-sw_heating_var_n;

%% make plot

figure('position',[100 100 700 300]); 
subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
plot(pw(:,1),fluxes_var_rs(:,3),'b-','Linewidth',1);
plot(pw(:,1),fluxes_var_rs(:,2),'r-','Linewidth',1);
plot(pw(:,1),QL_fixn,'k--','Linewidth',1);
plot(pw(:,1),Q_fixn,'k-','Linewidth',1);
[~,ind_max]=max(QL_fixn);
plot(pw(ind_max,1),QL_fixn(ind_max),'+','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixn);
plot(pw(ind_max,1),Q_fixn(ind_max),'+','Linewidth',1,'color',[0 0 0]);
xlim([0 90]);
ylim([0 500]);
box on;
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
ylabel('Radiative flux [W m$^{-2}$]','Interpreter','Latex');
LEG=legend('SLW','OLR','$\hat{Q}_L$','$\hat{Q} = \hat{Q}_L+\hat{Q}_S$','Location','NorthEast');
set(LEG,'Interpreter','latex')
title('a) Varied $r_s$','Fontweight','normal','Interpreter','Latex');
set(LEG,'position', [0.24 0.62 0.22 0.27]);

subplot(1,2,2);
hold on;
plot(pw(:,2),fluxes_var_n(:,3),'b-','Linewidth',1);
plot(pw(:,2),fluxes_var_n(:,2),'r-','Linewidth',1);
plot(pw(:,2),QL_fixrs,'k--','Linewidth',1);
plot(pw(:,2),Q_fixrs,'k-','Linewidth',1);
[~,ind_max]=max(QL_fixrs);
plot(pw(ind_max,1),QL_fixrs(ind_max),'+','Linewidth',1,'color',[0 0 0]);
[~,ind_max]=max(Q_fixrs);
plot(pw(ind_max,1),Q_fixrs(ind_max),'+','Linewidth',1,'color',[0 0 0]);
xlim([0 90]);
ylim([0 500]);
box on;
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
%ylabel('Radiative flux (W m$^{-2}$)','Interpreter','Latex');
%legend('OLR','SLF','Q_{LW}','Q_{LW+SW}','Location','best');
title('b) Varied $n$','Fontweight','normal','Interpreter','Latex');

gcfsavepdf([basedir 'Figure08_LW_SW_rrtmg.pdf']);
 
