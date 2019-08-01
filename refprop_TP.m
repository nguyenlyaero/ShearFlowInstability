function prop = refprop_TP(id, T, P, Y, paramArray)
% Calculate non-dim thermo property from non-dim T and P
    [T_base, P_ref, rho_ref, Rgas, cp_ref, mu_ref, lambda_ref, Ma, Pr, Rey, h_ref, T_ref, d_ref, Pe] = deal(paramArray{:});
    
    M_CH4 = 16.043; % kg/kmol
    M_O2 = 32.00; % kg/mol
    x_CH4 = M_O2*Y/(M_CH4*(1-Y)+M_O2*Y);
    M = x_CH4*M_CH4 + (1-x_CH4)*M_O2;
    %str = sprintf('CH4[%.12f]&O2[%.12f]',x_CH4,1-x_CH4);
    switch 1
        case strcmp(id, 'rho')
            %prop = PropsSI('D', 'T', T*T_ref, 'P', P*P_ref, str)/rho_ref;
            rho_hat = PC_SAFT_PT_Cubic(P*P_ref, T*T_ref ,x_CH4, 'CH4', 'O2', 0.5); % kmol/m^3
            prop = rho_hat*M/rho_ref;
    end
end