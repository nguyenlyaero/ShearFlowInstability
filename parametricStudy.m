function alphaMax = parametricStudy(P_Pc, Rey, y_pb, alphaLeft, alphaRight)
    %% Init
    N = 800;
    Rbar = 8.3144598; % J/mol-K

    P_ref = P_Pc*50.43e5; % Pa

    M_CH4 = 16.043; % kg/kmol
    M_O2 = 32.00; % kg/mol
    x_CH4_ref = M_O2*0.5/(M_CH4*0.5+M_O2*0.5);

    M = x_CH4_ref*M_CH4 + (1-x_CH4_ref)*M_O2;
    Rgas = Rbar/M*1000; % J/kg-K
    
    % Pseudoboiling stuffs
    Y_CH4_pb = 0.5*(1 + tanh(y_pb));
    x_CH4_pb = M_O2*Y_CH4_pb/(M_CH4*(1-Y_CH4_pb)+M_O2*Y_CH4_pb);
    T_pb = fminsearch(@(T) ObjectiveFun(T, P_ref, x_CH4_pb), 200);
    
    T_base = T_pb/(1 + 0.2*tanh(y_pb));
    
    % Derived thermodynamic
    rho_ref = PC_SAFT_PT_Cubic(P_ref,T_base,x_CH4_ref,'CH4','O2',0.5)*M; % kg/m^3
    [~, cp_ref] = IdealEnthalpy(T_base, x_CH4_ref, 'CH4', 'O2'); 
    cp_ref = cp_ref*Rgas; % J/kg-K
    h_ref = (cp_ref/Rgas*P_ref/rho_ref); % J/kg

    [mu_ref, lambda_ref] = Chung_PT(P_ref, T_base, x_CH4_ref, 'CH4', 'O2'); % Pa-s % W/m-k
    d_ref = Diffusivity_PT(P_ref, T_base, x_CH4_ref, 'CH4', 'O2'); % m^2/s

    T_ref = P_ref/Rgas/rho_ref; % K

    Pr = mu_ref*cp_ref/lambda_ref;
    
    % Nondimentional stuffs
    height = 0.3e-3; %m
    U_ref = Rey*(mu_ref/height/rho_ref);


    Ma = U_ref/sqrt(P_ref/rho_ref);

    Pe = U_ref*height/d_ref;

    paramArray = {T_base P_ref rho_ref Rgas cp_ref mu_ref lambda_ref Ma Pr Rey h_ref T_ref d_ref Pe};
    
    % Generate base flow
    [u0__x_vec, du0x_dxi_vec, d2u0x_dxidxi_vec,...
        rho0_vec, drho0_dxi_vec, d2rho0_dxidxi_vec,...
        T0_vec, dT0_dxi_vec, d2T0_dxidxi_vec, ...
        Y0_vec, dY0_dxi_vec, d2Y0_dxidxi_vec, ...
        D12_vec, dD_dT_vec, dD_drho_vec, dD_dY_vec, d2D_dTdT_vec, d2D_dTdrho_vec, d2D_dTdY_vec, d2D_drhodrho_vec, d2D_drhodY_vec, d2D_dYdY_vec,...
        mu_vec, dmu_dT_vec, dmu_drho_vec, dmu_dY_vec, d2mu_dTdT_vec, d2mu_dTdrho_vec, d2mu_dTdY_vec, d2mu_drhodrho_vec, d2mu_drhodY_vec, d2mu_dYdY_vec, ...
        P_vec, dp_dT_vec, dp_drho_vec, dp_dY_vec, d2p_dTdT_vec, d2p_dTdrho_vec, d2p_dTdY_vec, d2p_drhodrho_vec, d2p_drhodY_vec, d2p_dYdY_vec, ...
        h_vec, dh_dT_vec, dh_drho_vec, dh_dY_vec,...
        lambda_vec, dlambda_dT_vec, dlambda_drho_vec, dlambda_dY_vec, d2lambda_dTdT_vec, d2lambda_dTdrho_vec, d2lambda_dTdY_vec, d2lambda_drhodrho_vec, d2lambda_drhodY_vec, d2lambda_dYdY_vec] = Calc_BaseFlow_Nodes(N, paramArray);
    
    % Solve for alphaMax
    alphaMax = getAlphaMax(N, paramArray, alphaLeft, alphaRight, ...
    u0__x_vec, du0x_dxi_vec, d2u0x_dxidxi_vec,...
        rho0_vec, drho0_dxi_vec, d2rho0_dxidxi_vec,...
        T0_vec, dT0_dxi_vec, d2T0_dxidxi_vec, ...
        Y0_vec, dY0_dxi_vec, d2Y0_dxidxi_vec, ...
        D12_vec, dD_dT_vec, dD_drho_vec, dD_dY_vec, d2D_dTdT_vec, d2D_dTdrho_vec, d2D_dTdY_vec, d2D_drhodrho_vec, d2D_drhodY_vec, d2D_dYdY_vec,...
        mu_vec, dmu_dT_vec, dmu_drho_vec, dmu_dY_vec, d2mu_dTdT_vec, d2mu_dTdrho_vec, d2mu_dTdY_vec, d2mu_drhodrho_vec, d2mu_drhodY_vec, d2mu_dYdY_vec, ...
        P_vec, dp_dT_vec, dp_drho_vec, dp_dY_vec, d2p_dTdT_vec, d2p_dTdrho_vec, d2p_dTdY_vec, d2p_drhodrho_vec, d2p_drhodY_vec, d2p_dYdY_vec, ...
        h_vec, dh_dT_vec, dh_drho_vec, dh_dY_vec,...
        lambda_vec, dlambda_dT_vec, dlambda_drho_vec, dlambda_dY_vec, d2lambda_dTdT_vec, d2lambda_dTdrho_vec, d2lambda_dTdY_vec, d2lambda_drhodrho_vec, d2lambda_drhodY_vec, d2lambda_dYdY_vec);
end

