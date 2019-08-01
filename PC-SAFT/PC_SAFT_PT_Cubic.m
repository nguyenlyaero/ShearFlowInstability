function [rho_hat,ah, log_phi_vec, h_hat] = PC_SAFT_PT_Cubic(P,T,x,S1,S2,eta0)
% a: reduced Helholtz free energy (a/RT)
% rho_hat: molar density (kmol/m^3)
%
% T: [K]
% P: Pressure (Pa)
% x: molar fraction of first species
% str: species (e.g. 'H2&O2')
% eta0: initial guess, 1e-10 for vapor, 0.5 for liquid

    MaxIter = 1000;
    RelTol = 1e-10;

    N_Av=6.022e23;
    %% Data
    %% Get Species Data
    if nargin == 5
        [Tc, Pc, acentric, M, m , epsilon_k, sigma] = getSpecies(S1);
        m_vec = m;
        epsilon_k_vec = epsilon_k;
        sigma_vec = sigma;
        x_vec = 1;
        k_bin = 0;
        
        eta0 = S2;
        
        flag = 1;
    else % nargin == 6
        [Tc1, Pc1, acentric1, M1, m1 , epsilon_k1, sigma1] = getSpecies(S1);
        [Tc2, Pc2, acentric2, M2, m2 , epsilon_k2, sigma2] = getSpecies(S2);
        m_vec = [m1; m2];
        epsilon_k_vec = [epsilon_k1; epsilon_k2];
        sigma_vec = [sigma1; sigma2];
        x_vec = [x; 1-x];
        k_bin = getBinaryInteraction(S1, S2);
        
        flag = 0;
    end
    
    Nspecies = length(m_vec);

    %%
    m_bar = sum(x_vec.*m_vec);
    d_vec = sigma_vec.*(1-0.12.*exp(-3*epsilon_k_vec/T)); % Temp. Dependent segment diameter
    
    %% Solve for eta
    
    function [Perr,ah,log_phi_vec,h_hat] = objective_fun(eta)
        if flag == 1
            step = 1e-8;
            rho_hati = eta*(pi/6*sum(x_vec.*m_vec.*d_vec.^3))^(-1)/N_Av*1e30*1e-3;
            [Pi,ah,log_phi_vec,h_hat] = PC_SAFT_DT_Cubic(rho_hati,T,x,S1);
            Perr = (Pi-P)/P;
        else % nargin == 6
            step = 1e-8;
            rho_hati = eta*(pi/6*sum(x_vec.*m_vec.*d_vec.^3))^(-1)/N_Av*1e30*1e-3;
            [Pi,ah,log_phi_vec,h_hat] = PC_SAFT_DT_Cubic(rho_hati,T,x,S1,S2);
            Perr = (Pi-P)/P;
        end
    end

    function dPerr_deta = objective_derivative(eta)
        step = 1e-30;
        dPerr_deta = imag(objective_fun(eta*complex(1,step)))/eta/step;
    end
    
    eta = eta0;
    [Perr,ah,log_phi_vec,h_hat] = objective_fun(eta);
    dPerr_deta = objective_derivative(eta);
    iter = 0;
    
    while iter < MaxIter && abs(Perr) > RelTol
        eta = eta - Perr/dPerr_deta;
        [Perr,ah,log_phi_vec,h_hat] = objective_fun(eta);
        dPerr_deta = objective_derivative(eta);
        iter = iter+1;
    end
        
    assert(abs(Perr) < RelTol);
    
    rho_hat = eta*(pi/6*sum(x_vec.*m_vec.*d_vec.^3))^(-1)/N_Av*1e30*1e-3;
end

