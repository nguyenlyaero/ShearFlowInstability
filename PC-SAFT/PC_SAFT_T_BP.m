function [P, y] = PC_SAFT_T_BP(T,x,S1,S2,P0,y0)
% a: reduced Helholtz free energy (a/RT)
% rho_hat: molar density (kmol/m^3)
%
% T: [K]
% P: Pressure (Pa)
% z: molar fraction of first species
% str: species (e.g. 'H2&O2')
% eta0: initial guess, 1e-10 for vapor, 0.5 for liquid

    MaxIter = 1000;
    RelTol = 1e-8;

    N_Av=6.022e23;
    %% Get Species Data
    [Tc1, Pc1, acentric1, M1, m1 , epsilon_k1, sigma1] = getSpecies(S1);
    [Tc2, Pc2, acentric2, M2, m2 , epsilon_k2, sigma2] = getSpecies(S2);
    m_vec = [m1; m2];
    epsilon_k_vec = [epsilon_k1; epsilon_k2];
    sigma_vec = [sigma1; sigma2];
    x_vec = [x; 1-x];
    k_bin = getBinaryInteraction(S1, S2);
    
    Pc_vec = [Pc1; Pc2];
    Tc_vec = [Tc1; Tc2];
    omega_vec = [acentric1; acentric2];

    Nspecies = length(m_vec);
    

    %%
    % Initial Guess
    if (~exist('y0','var') || ~exist('P0','var'))
        P_sat_vec = Pc_vec.*10.^(7/3*(1+omega_vec).*(1-Tc_vec/T));
        P0 = sum(x_vec.*P_sat_vec);
        y0_vec = x_vec(1)*P_sat_vec(1)/P0;
        y0_vec = [y0_vec; 1-y0_vec];
    else
        y0_vec = [y0; 1-y0];
    end
    % first iteration    
    P = P0;
    y_vec = y0_vec;
    
    [~,~,log_phi_vec_VAPOR] = PC_SAFT_PT_Cubic(P,T,y_vec(1),S1,S2,1e-10);
    phi_vec_VAPOR = exp(log_phi_vec_VAPOR);

    [~,~,log_phi_vec_LIQUID] = PC_SAFT_PT_Cubic(P,T,x_vec(1),S1,S2,0.5);
    phi_vec_LIQUID = exp(log_phi_vec_LIQUID);   
    
    K_vec = phi_vec_LIQUID./phi_vec_VAPOR;
    yT = sum(x_vec.*K_vec);
    
    yT_oldold = yT;
    P_oldold = P;
    
    % second iteration
    P = P*yT;
    y_vec = x_vec(1)*K_vec(1)/yT;
    y_vec = [y_vec; 1-y_vec];

    [~,~,log_phi_vec_VAPOR] = PC_SAFT_PT_Cubic(P,T,y_vec(1),S1,S2,1e-10);
    phi_vec_VAPOR = exp(log_phi_vec_VAPOR);

    [~,~,log_phi_vec_LIQUID] = PC_SAFT_PT_Cubic(P,T,x_vec(1),S1,S2,0.5);
    phi_vec_LIQUID = exp(log_phi_vec_LIQUID);   
    
    K_vec = phi_vec_LIQUID./phi_vec_VAPOR;
    yT = sum(x_vec.*K_vec);

    yT_old = yT;
    P_old = P;
    
    % if we need third and more iterations
    while (abs(yT-1)>1e-10)

        P = P_oldold + (1-yT_oldold)/(yT_old-yT_oldold)*(P_old-P_oldold);
        y_vec = x_vec(1)*K_vec(1)/yT;
        y_vec = [y_vec; 1-y_vec];
        
        [~,~,log_phi_vec_VAPOR] = PC_SAFT_PT_Cubic(P,T,y_vec(1),S1,S2,1e-10);
        phi_vec_VAPOR = exp(log_phi_vec_VAPOR);

        [~,~,log_phi_vec_LIQUID] = PC_SAFT_PT_Cubic(P,T,x_vec(1),S1,S2,0.5);
        phi_vec_LIQUID = exp(log_phi_vec_LIQUID);   

        K_vec = phi_vec_LIQUID./phi_vec_VAPOR;
        yT = sum(x_vec.*K_vec);

        yT_oldold = yT_old;
        P_oldold = P_old;
        yT_old = yT;
        P_old = P;
    end
    y = y_vec(1);
end

