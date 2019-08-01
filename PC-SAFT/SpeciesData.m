% Species data
% Senol 2011
Tc_H2 = 33.2; % K
Pc_H2 = 13.0e5; % Pa
acentric_H2 = -0.218;
M_H2 = 2; % [g/mol]
m_H2 = 0.935864;
epsilon_k_H2 = 25.62934; % [K]
sigma_H2 = 2.912599; % [Angstrom]

% Stoll Vrabec Hasse 2003
Tc_O2 = 154.6; % K
Pc_O2 = 50.4e5; % Pa
acentric_O2 = 0.025;
M_O2 = 32; % [g/mol]
m_O2 = 1.1457;
epsilon_k_O2 = 113.98; % [K]
sigma_O2 = 3.1711; % [Angstrom]

% Gross Sadowski 2001
Tc_CH4 = 190.6; % K
Pc_CH4 = 46.0e5; % Pa
acentric_CH4 = 0.008;
M_CH4 = 16.043; % [g/mol]
m_CH4 = 1; % # Segments per chain
sigma_CH4 = 3.7039; % Segment Diameter [Angstrom]
epsilon_k_CH4 = 150.03; % Potential Depth [K]

Tc_CO2 = 304.2; % K
Pc_CO2 = 73.8e5; % Pa
acentric_CO2 = 0.224;
M_CO2=44.01; % [g/mol]
m_CO2=2.0729; % # Segments per chain
sigma_CO2=2.7852; % Segment Diameter [Angstrom]
epsilon_k_CO2=169.21; % Potential Depth [K]

% Delters 1993
k_H2_O2 = 1-0.92683;

save('SpeciesData.mat');

% Lease Squares with CoolProp
k_CH4_O2 = 0.098561433743562;