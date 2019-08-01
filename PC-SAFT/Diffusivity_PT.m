function [Diffusivity] = Diffusivity_PT(P,T,x,S1,S2)
% Calculate molecular diffusivity of binary mixture using Chapman-Enskog
% molecular dynamics model with Takahashi's high pressure correction
% Input: Pressure (Pa)

k_b = 1.380649e-23; % Boltzmann Constant (J/K)
N_Av=6.022e23;
% Get species data
[Tc1, Pc1, ~, M1, ~ , epsilon_k1, ~, sigma_coll1] = getSpecies(S1);
[Tc2, Pc2, ~, M2, ~ , epsilon_k2, ~, sigma_coll2] = getSpecies(S2);

m_1 = M1/N_Av;
m_2 = M2/N_Av;
sigma_1 = sigma_coll1*1e-10; % [m]
sigma_2 = sigma_coll2*1e-10; % [m]

% Mixing Rule:
m_12 = m_1*m_2/(m_1+m_2); % Molecular mass (kg)
sigma_12 = (sigma_1 + sigma_2)/2; % Collision diameter [m]
epsilon_k12 = sqrt(epsilon_k1*epsilon_k2); % Interaction well depth (J)

Pcm = x*Pc1 + (1-x)*Pc2;
Tcm = x*Tc1 + (1-x)*Tc2;

% Chapman-Enskog
a1 = 1.0548; a2 = 0.1550; a3 = 0.55909; a4 = 2.1705;
T_star = T/epsilon_k12; % Reduce temperature
Omega_12 = a1*T_star^(-a2) + (T_star + a3)^(-a4); % Collision Integral, Lennard-Jones 12-6 Potential

D_l = 3/16*sqrt(2*pi*k_b^3*T^3/m_12)/(P*pi*sigma_12^2*Omega_12);

Diffusivity = D_l;

% Takahashi correction
% P_r_vec = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0 4.0 5.0];
% DP_vec = [1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.02 1.02 1.02 1.02 1.03 1.03 1.04 1.05 1.06 1.07];
% A_vec = [0.038042 0.067433 0.098317 0.137610 0.175081 0.216376 0.314051 0.385736 0.514553 0.599184 0.557725 0.593007 0.696001 0.790770 0.502100 0.837452 0.890390];
% B_vec = [1.52267 2.16794 2.42910 2.77605 2.98256 3.11384 3.50264 3.07773 3.54744 3.61216 3.41882 3.18415 3.37660 3.27984 2.39031 3.23513 3.13001];
% C_vec = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.141211 0.278407 0.372683 0.504894 0.678469 0.665702 0.0 0.602907 0.0 0.0];
% E_vec = [13.45454 13.45454 13.45454 13.45454 13.45454 13.45454 13.45454 13.45454 14.00000 10.00900 8.57519 10.37483 11.21674 11.21674 6.19043 6.19043 6.19043];
% 
% P_r = P/Pcm;
% T_r = T/Tcm;
% 
% assert(P_r<5);
% if P_r<0.1
%     Diffusivity = D_l;
% else
%     for idx=1:length(P_r_vec)
%         if P_r_vec(idx)>P_r
%             frac = (P_r-P_r_vec(idx-1))/(P_r_vec(idx)-P_r_vec(idx-1));
%             break;
%         end
%     end
% 
%     % Corr1    
%     DP = DP_vec(idx);
%     A = A_vec(idx);
%     B = B_vec(idx);
%     C = C_vec(idx);
%     E = E_vec(idx);
%     f = C*T_r^(-E);
%     if f>1 
%         f = 1;
%     end
% 
%     Corr1 = DP*(1 - A*T_r^(-B))*(1-f);
%     
%     % Corr2    
%     DP = DP_vec(idx-1);
%     A = A_vec(idx-1);
%     B = B_vec(idx-1);
%     C = C_vec(idx-1);
%     E = E_vec(idx-1);
%     f = C*T_r^(-E);
%     if f>1 
%         f = 1;
%     end
% 
%     Corr2 = DP*(1 - A*T_r^(-B))*(1-f);
%     
%     % Interpolate
%     Corr = Corr2*(1-frac) + Corr1*frac;
%     Diffusivity = D_l*Corr;
% 
% end

end

