%% Figure 10

% Plots the mixing ratio, relative humidity, net upwards longwave and
% shortwave fluxes for Chuuk Lagoon's January 1996/1997/.../2016 soundings
% For each case, plot the moist cases (+1STD) in blue, 
% the average cases in green and the dry cases (-1STD) in red.

%% 1. Load data

load('Figure10_data.mat'); % Contains p[pressure used for radiative fluxes],
% pRH [pressure used for moisture variables], 
% mixing ratios r, relative humidities RH, 
% net upwards longwave flux LW and net upwards shortwave flux SW.

%% 2. Plot

figure('position',[100 100 840 300]);
fsize = 10; %font size for all plots

% colors for dry, average, and moist soundings
cdry    = [1 0.5 0];
cav     = [0.5 0.5 0.5];
cmoist  = [0 0 0.75];

HP1=subplot(1,4,1);
hold on;
set(gca,'Fontsize',fsize);
plot(r_dry,pRH,'Linewidth',1,'color',cdry);
plot(r_av,pRH,'Linewidth',1,'color',cav);
plot(r_moist,pRH,'Linewidth',1,'color',cmoist);
box on;
set(gca,'Ydir','reverse');
ylabel('Pressure [hPa]','Interpreter','Latex');
xlabel('Mixing ratio, $r$ [kg/kg]','Interpreter','Latex');
title('a) Mixing ratio profile','FontWeight','normal','Interpreter','Latex');
L=legend('dry: $< -1\sigma$','mean profile','moist: $> +1\sigma$');

HP2=subplot(1,4,2);
hold on;
set(gca,'Fontsize',fsize);
plot(RH_dry,pRH,'Linewidth',1,'color',cdry);
plot(RH_av,pRH,'Linewidth',1,'color',cav);
plot(RH_moist,pRH,'Linewidth',1,'color',cmoist);
box on;
set(gca,'Ydir','reverse');
set(gca,'Yticklabel','');
xlabel('Relative humidity [\%]','Interpreter','Latex');
title('b) Relative humidity profile','FontWeight','normal','Interpreter','Latex');
%L=legend('<-1STD','Mean','>+1STD');
%set(L,'Location','NorthEast');

HP3=subplot(1,4,3);
hold on;
set(gca,'Fontsize',fsize);
plot(LW_dry,p,'Linewidth',1,'color',cdry);
plot(LW_av,p,'Linewidth',1,'color',cav);
plot(LW_moist,p,'Linewidth',1,'color',cmoist);
box on;
xlim([0 350]);
set(gca,'Ydir','reverse');
set(gca,'Yticklabel','');
xlabel('$(\mathcal{F}^{\uparrow}-\mathcal{F}^{\downarrow})_{\rm LW}$ [W m$^{-2}$]','Interpreter','Latex');
title('c) Net LW flux profile','FontWeight','normal','Interpreter','Latex');
%L=legend('<-1STD','Mean','>+1STD');
%set(L,'Location','NorthWest');

HP4=subplot(1,4,4);
hold on;
set(gca,'Fontsize',fsize);
plot(SW_dry,p,'Linewidth',1,'color',cdry);
plot(SW_moist,p,'Linewidth',1,'color',cmoist);
plot(SW_av,p,'Linewidth',1,'color',cav);
box on;
set(gca,'Ydir','reverse');
set(gca,'Yticklabel','');
xlabel('$(\mathcal{F}^{\downarrow}-\mathcal{F}^{\uparrow})_{\rm SW}$ [W m$^{-2}$]','Interpreter','Latex');
title('d) Net SW flux profile','FontWeight','normal','Interpreter','Latex');
%L=legend('<-1STD','Mean','>+1STD');
%set(L,'Location','NorthWest');

set(HP1,'position',[0.08  0.14  0.20  0.77]);
set(HP2,'position',[0.31  0.14  0.20  0.77]);
set(HP3,'position',[0.54  0.14  0.20  0.77]);
set(HP4,'position',[0.77  0.14  0.20  0.77]);
set(L,'Location','best','Interpreter','latex');

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure10_Chuuk.pdf']);