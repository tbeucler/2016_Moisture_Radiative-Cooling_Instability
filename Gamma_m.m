function [Gm] = Gamma_m(Ts,ps)

% Compute the average moist adiabatic lapse rate
% as a function of the surface temperature Ts [K]
% and the surface pressure ps [hPa]
% Accepts arrays if Ts and ps have the same length
% Follows Holton appendix D2
% Assumes hydrostasy and a constant "average moist adiabatic lapse rate"

%% 0. Constants

cp=1005; % Specific heat capacity of dry air at constant pressure
g=9.81; % Gravity constant
Gamma_d=g/cp; % Dry adiabatic lapse rate
Gm_range=linspace(Gamma_d/10,Gamma_d,10000); % Possible range of moist adiabatic lapse rate
Lv=2.5e6; % Latent heat of vaporization of water vapor
pref=500; % Reference pressure for average moist adiabatic lapse rate [hPa]
Rd=287; % Specific gas constant of dry air
Rv=461; % Specific gas constant of water vapor
eps=Rd/Rv; % Ratio of the molecular weight of water vapor to the molecular weight of dry air

%% 1. Fits the best moist adiabatic lapse rate based on the atmospheric properties at pref

L=length(Ts);
Gm=zeros(L,1);

for ind=1:L
    
    Tsa=Ts(ind);
    psa=ps(ind);
    
    Ttest=Tsa.*(pref/psa).^(Gm_range*Rd/g); % Temperature at the reference pressure for different lapse rates
    rsat=r_sat(pref,Ttest); % Saturation specific humidity at the reference temp/pres for dif. lapse rates
    
    term_0=Gm_range./Gamma_d;
    term_1=1+Lv*rsat./(Rd*Ttest);
    term_2=1+eps*Lv^2.*rsat./(cp*Rd*Ttest.^2);
    
    [~,ind_min]=min(abs(term_1./term_2-term_0));
    Gm(ind)=Gm_range(ind_min);
    
end

end