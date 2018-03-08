%% Figure09_LW_SW_rrtmg_clouds
% script to plot tropospheric radiative cooling against column water vapor
% for rrtm calculations, with extremely idealized water vapor profiles
%     r(p) = r_s (p/p_s)^n
% column water vapor is controlled by varying either r_s or n
%
% this script looks across several different cloud profile assumptions to
% see how clouds alter importance of vertical profile of water vapor
% perturbation

thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
datadir = [basedir(1:strfind(basedir,'Figure_scripts')-1) 'rrtmg_idealized/'];

% load precipitable water arrays for two variants (var-rs and var-n)
load([datadir 'pw_idealized_chuuk.mat']);

cloud_type = {''; 'low-thin-'; 'low-thick-'; 'high-thin-'; 'high-thick-'};

% array of total atm. cooling as function of PW, moisture perturbation
% type, and cloud type
qtot_atm = zeros(50,2,5);

% scale factor for sw fluxes (should depend on zenith angle used and latitude)
sw_scale_factor = 4/pi^2;
% 4/pi^2 appropriate for equatorial equinox with cos(zeta)=pi/4

% loop over input files/cloud type cases
for ic=1:length(cloud_type)
    fluxes_var_rs = load([datadir 'output_summary_lw-sw_' char(cloud_type(ic)) 'var_rs.txt']);
    fluxes_var_n  = load([datadir 'output_summary_lw-sw_' char(cloud_type(ic)) 'var_n.txt']);

    qtot_atm(:,1,ic) = fluxes_var_rs(:,2)-fluxes_var_rs(:,3) -...
        sw_scale_factor*(fluxes_var_rs(:,4)-fluxes_var_rs(:,5));
    
    qtot_atm(:,2,ic) = fluxes_var_n(:,2)-fluxes_var_n(:,3) -...
        sw_scale_factor*(fluxes_var_n(:,4)-fluxes_var_n(:,5));
end

%% make plot

cmap = lines(5);

figure('position',[100 100 700 300]); 
subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
for icloud = 1:5
    plot(pw(:,1),qtot_atm(:,1,icloud),'Linewidth',1);
end
for icloud = 1:5
    [~,ind_max]=max(qtot_atm(:,1,icloud));
    plot(pw(ind_max,1),qtot_atm(ind_max,1,icloud),'+','Linewidth',1,'color',cmap(icloud,:));
end
xlim([0 90]);
ylim([0 400]);
box on;
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
ylabel('Radiative cooling, $\hat{Q}$ [W m$^{-2}$]','Interpreter','Latex');
LEG=legend('No clouds','Low thin clouds','Low thick clouds','High thin clouds','High thick clouds','Location','best');
title('a) Varied $r_s$','Fontweight','normal','Interpreter','latex');

subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
for icloud = 1:5
    plot(pw(:,2),qtot_atm(:,2,icloud),'Linewidth',1);
end
for icloud = 1:5
    [~,ind_max]=max(qtot_atm(:,2,icloud));
    plot(pw(ind_max,2),qtot_atm(ind_max,2,icloud),'+','Linewidth',1,'color',cmap(icloud,:));
end
xlim([0 90]);
ylim([0 400]);
box on;
xlabel('Column water vapor, $\hat{r}$ [kg m$^{-2}$]','Interpreter','Latex');
%ylabel('Atmospheric LW+SW cooling (W m$^{-2}$)','Interpreter','Latex');
%legend('No clouds','Low thin clouds','Low thick clouds','High thin clouds','High thick clouds','Location','best');
title('b) Varied $n$','Fontweight','normal','Interpreter','Latex');

gcfsavepdf([basedir 'Figure09_LW_SW_rrtmg_clouds.pdf']);
 
