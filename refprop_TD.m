function prop = refprop_TD(id, T, rho, Y, paramArray)
% Calculate non-dim thermo property from non-dim T and rho
    [T_base, P_ref, rho_ref, Rgas, cp_ref, mu_ref, lambda_ref, Ma, Pr, Rey, h_ref, T_ref, d_ref, Pe] = deal(paramArray{:});
    
    M_CH4 = 16.043; % kg/kmol
    M_O2 = 32.00; % kg/mol
    x_CH4 = M_O2*Y/(M_CH4*(1-Y)+M_O2*Y);
    M = x_CH4*M_CH4 + (1-x_CH4)*M_O2; % kg/kmol
    %str = sprintf('CH4[%.12f]&O2[%.12f]',x_CH4,1-x_CH4);
    switch 1
        case strcmp(id, 'mu')
            rho_hat = rho*rho_ref/M;
            P_dim = PC_SAFT_DT_Cubic(rho_hat, T*T_ref, x_CH4, 'CH4', 'O2');
            prop = Chung_PT(P_dim, T*T_ref, x_CH4, 'CH4', 'O2')/mu_ref;
        case strcmp(id, 'P')
            rho_hat = rho*rho_ref/M;
            prop = PC_SAFT_DT_Cubic(rho_hat, T*T_ref, x_CH4, 'CH4', 'O2')/P_ref;
        case strcmp(id, 'h')
            rho_hat = rho*rho_ref/M;
            [~,~,~,h_hat] = PC_SAFT_DT_Cubic(rho_hat, T*T_ref, x_CH4, 'CH4', 'O2');
            prop = h_hat/M*1000/h_ref;
        case strcmp(id, 'lambda')
            rho_hat = rho*rho_ref/M;
            P_dim = PC_SAFT_DT_Cubic(rho_hat, T*T_ref, x_CH4, 'CH4', 'O2');
            [~, prop] = Chung_PT(P_dim, T*T_ref, x_CH4, 'CH4', 'O2');
            prop = prop/lambda_ref;
        case strcmp(id, 'd')
            rho_hat = rho*rho_ref/M;
            P_dim = PC_SAFT_DT_Cubic(rho_hat, T*T_ref, x_CH4, 'CH4', 'O2');
            prop = Diffusivity_PT(P_dim, T*T_ref, x_CH4, 'CH4', 'O2')/d_ref;           
    end
end