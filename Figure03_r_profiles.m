%% Figure 3

% Mixing ratio profile = r(p)=rs*(p/ps)^n perturbed in two ways
% Varying rs while n is fixed (varying surface mixing ratio, fixed shape)
% Varying n while rs is fixed (varying shape, fixed surface mixing
% ratio)
% We use constant increments of column water vapor:
% CWV = int(r(p)dp) = rs*ps/((n+1)*g)

p           =   linspace(0,1000,1000);  % Pressure
Lp          =   length(p);

n           =   3;                      % Best-fit exponent value from Chuuk Lagoon
profile_r   =   (p/1000).^n;            % Mixing ratio profile

Lvary       =   5;                      % Number of profiles

rsa         =   linspace(0.5,1.5,Lvary);% Varying surface mixing ratio
np1a        =   (n+1)./rsa;             % Equivalent increments of varying exponent n

profile_r_fixn=zeros(Lp,Lvary);
profile_r_fixrs=zeros(Lp,Lvary);

for irs=1:Lvary
    profile_r_fixn(:,irs)=rsa(irs).*(p/1000).^2; % Variation at fixed n
    profile_r_fixrs(:,irs)=(p/1000).^(np1a(irs)-1); % Variation at fixed r_s
end

% define colormap for plots
cmap    = [1 0.5 0; 0.75 0.5 0.25; 0.5 0.5 0.5; 0.25 0.25 1; 0 0 0.75];
%cmap    = flipud(parula(Lvary+1));
%cmap    = cmap(2:Lvary+1,:);

figure('position',[100 100 700 300]); 

subplot(1,2,1);
hold on;
set(gca,'Fontsize',11);
h1 = plot(profile_r_fixn,p,'Linewidth',1.5);
set(h1, {'color'}, num2cell(cmap, 2));
box on;
ylim([0 1000]);
xlim([0 1.5]);
set(gca,'Ydir','reverse')
xlabel('Normalized mixing ratio, $r(p)/r_{s0}$','Interpreter','Latex');
ylabel('Pressure [hPa]','Interpreter','Latex');
title('a) Varied $r_s$','Fontweight','normal','Interpreter','Latex');

subplot(1,2,2);
hold on;
set(gca,'Fontsize',11);
h2=plot(profile_r_fixrs,p,'Linewidth',1.5);
set(h2, {'color'}, num2cell(cmap, 2));
box on;
ylim([0 1000]);
xlim([0 1.5]);
set(gca,'Ydir','reverse')
xlabel('Normalized mixing ratio, $r(p)/r_{s0}$','Interpreter','Latex');
title('b) Varied $n$','Fontweight','normal','Interpreter','Latex');

% Trick to get current directory on different machines
thisfile    = which(mfilename);
basedir     = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'Figure03_r_profiles.pdf']);