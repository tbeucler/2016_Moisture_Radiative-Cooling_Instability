%% Figure 11

% Scatter plot of the longwave and total clear-sky MRCI growth rates
% Antilles (ANT) stations are marked by red dots, West Pacific (PAC)
% stations by green dots and Amazonian (AMA) stations by blue dots

%% 1. Load data

load('Figure11_data.mat'); % r is column water vapor and s is slope
% sLW is the slope when only the longwave feedback is taken into account
% sSW is the slope when only the shortwave feedback is taken into account

%% 2. Plot

figure('position',[100 100 700 300]);
marker_linewidth = 1;
marker_size = 5;

subplot(1,2,1)
hold on;
set(gca,'Fontsize',11);
plot(r_AMA,sLW_AMA,'o','Markerfacecolor',[0 0.4 0],'Markeredgecolor',[1 1 0],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
plot(r_ANT,sLW_ANT,'s','Markerfacecolor',[0 0 1],'Markeredgecolor',[1 0 0],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
plot(r_PAC,sLW_PAC,'p','Markerfacecolor',[1 1 1],'Markeredgecolor',[0 0 0.5],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
ylim([0 2]);
box on;
xlabel('Column water vapor [kg m$^{-2}$]','Interpreter','latex');
ylabel('Growth rate [month$^{-1}$]','Interpreter','latex');
title('a) LW MRCI growth rate','Interpreter','latex');

subplot(1,2,2)
hold on;
plot(r_AMA,sLW_AMA+sSW_AMA,'o','Markerfacecolor',[0 0.4 0],'Markeredgecolor',[1 1 0],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
plot(r_ANT,sLW_ANT+sSW_ANT,'s','Markerfacecolor',[0 0 1],'Markeredgecolor',[1 0 0],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
plot(r_PAC,sLW_PAC+sSW_PAC,'p','Markerfacecolor',[1 1 1],'Markeredgecolor',[0 0 0.5],'LineWidth',marker_linewidth,'MarkerSize',marker_size);
ylim([0 2]);
box on;
xlabel('Column water vapor [kg m$^{-2}$]','Interpreter','latex');
title('b) LW+SW MRCI growth rate','Interpreter','latex');

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure11_Obs_MRCI_growthrate.pdf']);
