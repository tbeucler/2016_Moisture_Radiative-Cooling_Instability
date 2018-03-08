%% FigXX_rrtmg_idealized_snowball
% script to plot tropospheric radiative cooling against column water vapor
% for rrtm calculations, with extremely idealized water vapor profiles
%     r(p) = r_s (p/p_s)^n
% column water vapor is controlled by varying either r_s or n
% AND very cold, snowball-earth-like soundings!
%

thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

% load precipitable water arrays for two variants (var-rs and var-n)
load([basedir 'rrtmg_idealized/pw_idealized_snowball.mat']);

% load fluxes for var-rs and var-n cases
fluxes_var_rs = load([basedir 'rrtmg_idealized/output_summary_snowball_lw-sw_var_rs.txt']);
fluxes_var_n = load([basedir 'rrtmg_idealized/output_summary_snowball_lw-sw_var_n.txt']);

% scale factor for sw fluxes (should depend on zenith angle used and latitude)
sw_scale_factor = 4/pi^2;
% 4/pi^2 appropriate for equatorial equinox with cos(zeta)=pi/4
sw_heating_var_rs = sw_scale_factor*(fluxes_var_rs(:,4)-fluxes_var_rs(:,5));
sw_heating_var_n = sw_scale_factor*(fluxes_var_n(:,4)-fluxes_var_n(:,5));

%% make plot

figure('position',[100 100 850 350]); 
subplot(1,2,1);
hold on;
plot(pw(:,1),fluxes_var_rs(:,2),'r-');
plot(pw(:,1),fluxes_var_rs(:,3),'-','Color',[1 0.5 0]);
plot(pw(:,1),fluxes_var_rs(:,2)-fluxes_var_rs(:,3),'b-');
plot(pw(:,1),fluxes_var_rs(:,2)-fluxes_var_rs(:,3)-sw_heating_var_rs,'b--');
xlim([0 2.5]);
ylim([0 250]);
box on;
xlabel('Column water vapor, $\hat{r}$ (kg m$^{-2}$)','Interpreter','Latex');
ylabel('Radiative flux (W m$^{-2}$)','Interpreter','Latex');
legend('OLR','SLF','Q_{LW}','Q_{LW+SW}','Location','best');
title('a) Varied r_s','Fontweight','normal');

subplot(1,2,2);
hold on;
plot(pw(:,2),fluxes_var_n(:,2),'r-');
plot(pw(:,2),fluxes_var_n(:,3),'-','Color',[1 0.5 0]);
plot(pw(:,2),fluxes_var_n(:,2)-fluxes_var_n(:,3),'b-');
plot(pw(:,2),fluxes_var_n(:,2)-fluxes_var_n(:,3)-sw_heating_var_n,'b--');

xlim([0 2.5]);
ylim([0 250]);
box on;
xlabel('Column water vapor, $\hat{r}$ (kg m$^{-2}$)','Interpreter','Latex');
ylabel('Radiative flux (W m$^{-2}$)','Interpreter','Latex');
legend('OLR','SLF','Q_{LW}','Q_{LW+SW}','Location','best');
title('b) Varied n_{\color{white} 0}','Fontweight','normal');

gcfsavepdf([basedir 'FigXX_rrtmg_lw-sw_idealized_snowball_var-rs-n.pdf']);
 
