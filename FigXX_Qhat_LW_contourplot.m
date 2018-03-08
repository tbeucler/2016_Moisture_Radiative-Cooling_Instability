%% FigXX_Qhat_LW_contourplot.m
% script to make a plot of total lw atmospheric cooling as a function of
% column-integrated water vapor (x-axis) and surface mixing ratio (y-axis)

% trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

% parameters/constants
%
% name      value       units           description
%--------------------------------------------------
k    =      0.1;        % m^2/kg        Water vapor window absorptiongray absorber coeff
TS   =      300;        % K             Surface temperature
Rdcp =      2/7;        % -             Ratio of Rd/Cp for use in calculating T-p profile
avm  =      0.6;        % -             Ratio of lapse rate to dry-adiabatic lapse rate
ps   =      10^5;       % Pa            Surface pressure
grav =      9.81;       % m/s^2         surface gravity
Dg   =      5/3;        % -             Gray radiation diffusivity factor

% temperature profile definition
pn   = (0:0.0002:1);
TA   = TS.*pn.^(Rdcp*avm);

num_rs      = 300;
num_rhat    = 300; 

% arrays of rs and rhat
rs   = 0.025*(1:num_rs)/num_rs;
rhat = (ps/grav)*max(rs)/2*(1:num_rhat)/num_rhat;

Qhat_atm    = zeros(length(rhat),length(rs));
n           = zeros(length(rhat),length(rs));

for i=1:length(rhat)
    for j=1:length(rs)       
        n(i,j) = (ps/grav)*rs(j)/rhat(i) - 1;
        if n(i,j)>=1
            tau0 = Dg*(ps*k/grav)*rs(j)/(n(i,j)+2);
            tau  = tau0*pn.^(n(i,j)+2);
            [lwnt, lwns] = OLR_band( tau, TA, 1);
            
            Qhat_atm(i,j) = lwnt-lwns;
        else
            Qhat_atm(i,j) = NaN;
        end
    end
end

%% make plot

% basic contouring of Qhat against rhat and rs
figure;
hc=imagesc(rhat,1000*rs,(Qhat_atm')/(5.67e-8*300^4));
set(hc,'alphadata',~isnan(Qhat_atm'));
set(gca,'Ydir','normal');
hold on;
colorbar;
axis square;

% add curve connecting maxima of Qhat at each value of rs
rhat_of_max_Qhat = zeros(length(rs),1);
for j=1:length(rs);
    max_Qhat            = max(Qhat_atm(:,j));
    rhat_of_max_Qhat(j) = rhat(Qhat_atm(:,j)==max_Qhat);
end
plot(rhat_of_max_Qhat,1000*rs,'r-');
ht = text(0.045*max(rhat),1000*0.025*max(rs),'Red line tracks max($\hat{Q}_{LW}(\hat{r})$) for each value of $r_s$','Color',[1 0 0]);
set(ht,'Interpreter','latex')

% add lines for constant n and label them
n_values    = [1 2 4 8 16];
num_n_lines = length(n_values);
rs_const_n    = zeros(length(rhat),num_n_lines);
for nn=1:num_n_lines;
    rs_const_n(:,nn) = rhat*(n_values(nn)+1)/(ps/grav);
    plot(rhat,1000*rs_const_n(:,nn),'Color','k');
    text_y      = 1000*0.92*max(rs);
    text_x      = 0.92*max(rs)*(ps/grav)/(n_values(nn)+1);
    text_angle  = atand((n_values(nn)+1)/2); 
    ht = text(text_x, text_y, sprintf('$n=%d$',n_values(nn)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Rotation',text_angle);
    set(ht,'Interpreter','latex');
end

% add axis labels 
xlabel('Column water vapor, $\hat{r}$ (kg/m$^2$)','Interpreter','latex');
ylabel('Surface mixing ratio, $r_s$ (g/kg)','Interpreter','latex');
title('Normalized column LW cooling, $\hat{Q}_{LW}/\sigma T_S^4$','Interpreter','latex','Fontweight','normal');
box on;
set(gca,'TickDir','out');

% save plot
gcfsavepdf([basedir 'FigXX_Qhat_LW_contourplot.pdf']);


