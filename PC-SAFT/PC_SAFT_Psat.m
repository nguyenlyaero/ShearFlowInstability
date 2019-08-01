function [Psat,rho_hat_v,rho_hat_l] = PC_SAFT_Psat(T,S1)
% PURE FLUIDS ONLY
% Psat: Saturation Pressure [Pa]
% rho_hat_v: vapor molar density [kmol/m^3]
% rho_hat_l: liquid molar density [kmol/m^3]
%
% T: [K]
% P: Pressure (Pa)
% str: species (e.g. 'H2', 'O2')
    MaxIter = 1000;
    RelTol = 1e-10;

    %% Get Species Data
    [Tc, Pc, acentric, M, m , epsilon_k, sigma] = getSpecies(S1);
    m_vec = m;
    epsilon_k_vec = epsilon_k;
    sigma_vec = sigma;
    x_vec = 1;
    k_bin = 0;
    
    Tr = T/Tc;
    
    %% Solve for Psat
    function [f_err, df_err_dP, rho_hat_v,rho_hat_l] = objective_fun(P)
        [rho_hat_v,~, log_phi_v] = PC_SAFT_PT_Cubic(P,T,1,S1,1e-10);
        [rho_hat_l,~, log_phi_l] = PC_SAFT_PT_Cubic(P,T,1,S1,0.5);
        f_err = (exp(log_phi_v)-exp(log_phi_l))/exp(log_phi_l);
        step = 1e-8;
        [~,~, log_phi_v_p] = PC_SAFT_PT_Cubic(P*(1+step),T,1,S1,1e-10);
        [~,~, log_phi_l_p] = PC_SAFT_PT_Cubic(P*(1+step),T,1,S1,0.5);
        f_err_p = (exp(log_phi_v_p)-exp(log_phi_l_p))/exp(log_phi_l_p);
        df_err_dP = (f_err_p-f_err)/(step*P);
    end

    Psat = Pc*10^(7/3*(1+acentric)*(1-1/Tr));
    [f_err, df_err_dP, rho_hat_v,rho_hat_l] = objective_fun(Psat);
    iter = 0;
    
    while iter < MaxIter && abs(f_err) > RelTol
        Psat = Psat - f_err/df_err_dP;
        [f_err, df_err_dP, rho_hat_v,rho_hat_l] = objective_fun(Psat);
        iter = iter+1;
    end
    
    assert(abs(f_err)<RelTol);
end


