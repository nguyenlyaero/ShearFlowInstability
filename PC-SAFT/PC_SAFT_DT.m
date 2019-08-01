function [P,a_h,log_phi_vec,h] = PC_SAFT_DT(rho_hat,T,x,S1,S2)
% a: reduced Helholtz free energy (a/RT)
% P: Pressure (Pa)
% log_phi_vec: log of fugacity coefficient for each species
%
% T: [K]
% rho_hat: molar density (kmol/m^3)
% x: molar fraction of first species
% str: species (e.g. 'H2&O2')

    N_Av=6.022e23;
    rho=rho_hat.*N_Av*10^3/(10^10)^3; % number molecular density (1/Angstrom)^3

    %% Get Species Data
    if nargin == 4
        [Tc, Pc, acentric, M, m , epsilon_k, sigma] = getSpecies(S1);
        m_vec = m;
        epsilon_k_vec = epsilon_k;
        sigma_vec = sigma;
        x_vec = 1;
        k_bin = 0;
    else % nargin == 5
        [Tc1, Pc1, acentric1, M1, m1 , epsilon_k1, sigma1] = getSpecies(S1);
        [Tc2, Pc2, acentric2, M2, m2 , epsilon_k2, sigma2] = getSpecies(S2);
        m_vec = [m1; m2];
        epsilon_k_vec = [epsilon_k1; epsilon_k2];
        sigma_vec = [sigma1; sigma2];
        x_vec = [x; 1-x];
        k_bin = getBinaryInteraction(S1, S2);
    end
    
    Nspecies = length(m_vec);
    
    %% Pressure
    if Nspecies == 1
        [Psat,rho_hat_v,rho_hat_l] = PC_SAFT_Psat(T,S1);
        if (rho_hat-rho_hat_v)*(rho_hat-rho_hat_l)<0 % VLE
            P = Psat;
        else
            [P,a_h,log_phi_vec,h] = PC_SAFT_DT_Cubic(rho_hat,T,x,S1);
        end
    else
        [P,a_h,log_phi_vec,h] = PC_SAFT_DT_Cubic(rho_hat,T,x,S1,S2);
    end
end

