%% FigXX_rrtmg_idealized
% script to plot tropospheric radiative cooling against column water vapor
% for rrtm calculations, with extremely idealized water vapor profiles
%     r(p) = r_s (p/p_s)^n
% column water vapor is controlled by varying either r_s or n

thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

% load precipitable water arrays for two variants (var-rs and var-n)
load([basedir 'rrtmg_idealized/pw_idealized_var-rs-n.mat']);

% load rrtmg_lw fluxes
lwfluxes_var_rs = load([basedir 'rrtmg_idealized/output_summary_lw_ts300-var_rs.txt']);
lwfluxes_var_n = load([basedir 'rrtmg_idealized/output_summary_lw_ts300-var_n.txt']);
% load rrtmg_sw fluxes
swfluxes_var_rs = load([basedir 'rrtmg_idealized/output_summary_sw_ts300-var_rs.txt']);
swfluxes_var_n = load([basedir 'rrtmg_idealized/output_summary_sw_ts300-var_n.txt']);

% scale factor for sw fluxes (should depend on zenith angle used and latitude)
sw_scale_factor = 3/8;
sw_heating_var_rs = sw_scale_factor*(swfluxes_var_rs(:,2)-swfluxes_var_rs(:,4));
sw_heating_var_n = sw_scale_factor*(swfluxes_var_n(:,2)-swfluxes_var_n(:,4));

%% make plot

figure('position',[100 100 850 350]); 
subplot(1,2,1);
hold on;
plot(pw(:,1),lwfluxes_var_rs(:,2),'r-');
plot(pw(:,1),lwfluxes_var_rs(:,4),'-','Color',[1 0.5 0]);
plot(pw(:,1),lwfluxes_var_rs(:,2)-lwfluxes_var_rs(:,4),'b-');
plot(pw(:,1),lwfluxes_var_rs(:,2)-lwfluxes_var_rs(:,4)-sw_heating_var_rs,'b--');
xlim([0 80]);
ylim([0 450]);
box on;
xlabel('Column water vapor, $\hat{r}$ (kg m$^{-2}$)','Interpreter','Latex');
ylabel('Radiative flux (W m$^{-2}$)','Interpreter','Latex');
legend('OLR','SLF','Q_{LW}','Q_{LW+SW}','Location','best');
title('a) Varied r_s','Fontweight','normal');

subplot(1,2,2);
hold on;
plot(pw(:,2),lwfluxes_var_n(:,2),'r-');
plot(pw(:,2),lwfluxes_var_n(:,4),'-','Color',[1 0.5 0]);
plot(pw(:,2),lwfluxes_var_n(:,2)-lwfluxes_var_n(:,4),'b-');
plot(pw(:,2),lwfluxes_var_n(:,2)-lwfluxes_var_n(:,4)-sw_heating_var_n,'b--');

xlim([0 80]);
ylim([0 450]);
box on;
xlabel('Column water vapor, $\hat{r}$ (kg m$^{-2}$)','Interpreter','Latex');
ylabel('Radiative flux (W m$^{-2}$)','Interpreter','Latex');
legend('OLR','SLF','Q_{LW}','Q_{LW+SW}','Location','best');
title('b) Varied n_{\color{white} 0}','Fontweight','normal');

gcfsavepdf([basedir 'FigXX_rrtmg_lw-sw_idealized_var-rs-n.pdf']);
 
