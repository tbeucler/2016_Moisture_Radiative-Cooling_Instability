%% Figure7

% Total (longwave-shortwave, clear-sky+cloudy)
% atmospheric radiative cooling as a function
% of the column water vapor at fixed surf. temp. Ts
% 1. At varying surface mixing ratio
% 2. At varying shape

%% 0. Constants and parameters

% 0.1 Constants

cp      =   1005;           % Specific heat capacity of dry air at constant pressure (J/kg/K)
D       =   1.66;           % Schwarzschild two-stream diffusivity factor
g       =   9.81;           % Surface gravity constant (m/s^2)
Gamma_d =   g/cp;           % Dry adiabatic lapse rate (K/m)
kLW     =   0.1;            % LW absorbtion coefficient of water vapor (kg^-1.m^2)
kSW     =   0.01;           % SW absorbtion coefficient of water vapor (kg^-1.m^2)
Lv      =   2.5e6;          % Latent heat of vaporization of water vapor (J/kg)
mu      =   pi/4;           % Insolation-weighted zenith angle at the Equator
n_RCE   =   3;              % Mixing ratio profile exponent based on Chuuk Lagoon
ps      =   1000;           % Surface pressure (hPa)
sig     =   5.67e-8;        % Stefan-Boltzmann constant (W/m^2/K^4)
S       =   sig*280^4;      % Insolation for an effective emission temperature of 280K (W/m^2)
Ts      =   300;            % Surface temperature (K)
Gm      =   Gamma_m(Ts,ps); % Average moist adiabatic lapse rate (K/m)
Rd      =   287;            % Dry specific gas constant (J/kg/K)
rs_RCE  =   r_sat(ps,Ts);   % Saturation mixing ratio 

% 0.2 Clear-sky parameters

alpha0      = 4*Gm*Rd/g;            % Exponent relating temperature and optical depth
alpha_RCE   = alpha0/(n_RCE+2);     % RCE value of this exponent
rh          = linspace(0.01,60,100);% Column water vapor
Lrh         = length(rh);
eps         = kSW/(kLW*mu*D);       % Ratio of shortwave to longwave optical depth

% 0.3 Cloud parameters

deltau      =   [0 0.5 5 0.5 5];                % Cloud's optical thickness
pc_ps       =   [0.900 0.900 0.900 0.200 0.200];% [Pressure level of the cloud]/[Surface pressure]
Lcloud      =   length(deltau);                 % Number of cloud situations

%% 1. Varying surface mixing ratio, fixed shape

tau_fixn=D*(n_RCE+1)*kLW/(n_RCE+2)*rh; % Optical depth of atm at fixed n
Q_fixn=zeros(Lrh,Lcloud); % Total atmospheric cooling

for irh=1:Lrh
    
    ta=tau_fixn(irh); % Clear-sky optical depth of the atmosphere
    
    for icloud=1:Lcloud
        
        pcps=pc_ps(icloud); % Ratio of the cloud pressure to the surface pressure
        del=deltau(icloud); % Optical thickness of the cloud
        tac=ta*pcps^(alpha0/alpha_RCE); % Level of the cloud in opt. th. coordinates
        
        y1_a = @(x) (x/ta).^alpha_RCE.*exp(x-ta-del);
        y2_a = @(x) (x/ta).^alpha_RCE.*exp(-x);
        y1_b = @(x) ((x-del)/ta).^alpha_RCE.*exp(x-ta-del);
        y2_b = @(x) ((x-del)/ta).^alpha_RCE.*exp(-x);
        Y1_a = integral(y1_a,0,tac);
        Y2_a = integral(y2_a,0,tac);
        Y1_b = integral(y1_b,tac+del,ta+del);
        Y2_b = integral(y2_b,tac+del,ta+del);
        
        Q_a=Y1_a+Y2_a; % Emission of the atmosphere above the cloud
        Q_b=Y1_b+Y2_b; % Emission of the atmosphere below the cloud
        Q_c=2*pcps^alpha0*sinh(del/2)*(exp(-tac-del/2)+exp(tac-ta-del/2)); % Emission of the cloud
        QL=-1+exp(-ta-del); % Clear-sky longwave surf/atm/space cooling
%        QS=(S/(sig*Ts^4))*(1-exp(-eps*(ta+del))); % Clear-sky shortwave heating
% twc - commented out, clouds absorb too much in SW in above formulation
        QS=(S/(sig*Ts^4))*(1-exp(-eps*ta)); % Clear-sky shortwave heating
        
        Q_fixn(irh,icloud)=sig*Ts^4*(QL-QS+Q_a+Q_b+Q_c); % Total atmospheric cooling
        
    end
    
end

%% 2. Varying shape, fixed surface mixing ratio

rsh=rs_RCE*100*ps/g; % Column water vapor if r(p)=rs
tau_fixrs=kLW*D*rsh*rh./(rsh+rh); % Optical depth of atm at fixed rs
alp_fixrs=alpha0*rh./(rsh+rh); % Exponent alpha at fixed rs
Q_fixrs=zeros(Lrh,Lcloud); % Total atmospheric cooling

for irh=1:Lrh
    
    ta=tau_fixrs(irh); % Clear-sky optical depth of the atmosphere
    al=alp_fixrs(irh); % Mixing ratio profile exponent
    
    for icloud=1:Lcloud
        
        pcps=pc_ps(icloud); % Ratio of the cloud pressure to the surface pressure
        del=deltau(icloud); % Optical thickness of the cloud
        tac=ta*pcps^(alpha0/al); % Level of the cloud in opt. th. coordinates
        
        y1_a = @(x) (x/ta).^al.*exp(x-ta-del);
        y2_a = @(x) (x/ta).^al.*exp(-x);
        y1_b = @(x) ((x-del)/ta).^al.*exp(x-ta-del);
        y2_b = @(x) ((x-del)/ta).^al.*exp(-x);
        Y1_a = integral(y1_a,0,tac);
        Y2_a = integral(y2_a,0,tac);
        Y1_b = integral(y1_b,tac+del,ta+del);
        Y2_b = integral(y2_b,tac+del,ta+del);
        
        Q_a=Y1_a+Y2_a; % Emission of the atmosphere above the cloud
        Q_b=Y1_b+Y2_b; % Emission of the atmosphere below the cloud
        Q_c=2*pcps^alpha0*sinh(del/2)*(exp(-tac-del/2)+exp(tac-ta-del/2)); % Emission of the cloud
        QL=-1+exp(-ta-del); % Clear-sky longwave surf/atm/space cooling
%        QS=(S/(sig*Ts^4))*(1-exp(-eps*(ta+del))); % Clear-sky shortwave heating
% twc - commented out, clouds absorb too much in SW in above formulation
        QS=(S/(sig*Ts^4))*(1-exp(-eps*ta)); % Clear-sky shortwave heating
       
        
        Q_fixrs(irh,icloud)=sig*Ts^4*(QL-QS+Q_a+Q_b+Q_c); % Total atmospheric cooling
        
    end
    
end

%% 3. Plot

figure('position',[100 100 700 300]);
cmap = lines(Lcloud);

subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
plot(rh,Q_fixn,'Linewidth',1);
[~,ind_max]=max(Q_fixn);
for icloud=1:Lcloud
    plot(rh(ind_max(icloud)),Q_fixn(ind_max(icloud),icloud),'+','Linewidth',1,'color',cmap(icloud,:));
end
box on;
ylim([-200 400]);
line([0 60],[0 0],'Color','k','LineStyle',:);
LEG=legend('No clouds','Low thin clouds','Low thick clouds','High thin clouds','High thick clouds');
xlabel('Column water vapor [kg m$^{-2}$]','Interpreter','latex');
ylabel('Radiative cooling, $\hat{Q}$ [W m$^{-2}$]','Interpreter','latex');
title('a) Varied $r_s$','Fontweight','normal','Interpreter','latex');
set(LEG,'Fontsize',10,'Location','Best');
set(LEG,'Position', [0.21 0.14 0.25 0.26]);
legend boxoff


subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
plot(rh,Q_fixrs,'Linewidth',1);
[~,ind_max]=max(Q_fixrs);
for icloud=1:Lcloud
    plot(rh(ind_max(icloud)),Q_fixrs(ind_max(icloud),icloud),'+','Linewidth',1,'color',cmap(icloud,:));
end
box on;
ylim([-200 400]);
line([0 60],[0 0],'Color','k','LineStyle',:);
%LEG=legend('No clouds','Low thin clouds','Low thick clouds','High thin clouds','High thick clouds');
xlabel('Column water vapor [kg m$^{-2}$]','Interpreter','latex');
title('b) Varied $n$','Fontweight','normal','Interpreter','latex');
%set(LEG,'Fontsize',10,'Location','Best');
%set(LEG,'Position', [0.61 0.165 0.25 0.25]);
%legend boxoff

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure07_LW_SW_gray_clouds.pdf']);