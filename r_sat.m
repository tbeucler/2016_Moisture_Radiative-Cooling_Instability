function r=r_sat(p,T)

% Computes the saturation mixing ratio [1] as a function of pressure [hPa]
% and temperature [K]
% Uses Bolton formula and Dalton's law

eps=0.621; % Ratio of the molecular weight of water vapor to the molecular weight of dry air 
x=T-273.15; % Converts temperature from [C] to [K]
e_sat=6.112.*exp(17.67.*x./(x+243.5)); % Bolton's formula
r=eps.*e_sat./(p-e_sat); % Dalton's law and definition of mixing ratio