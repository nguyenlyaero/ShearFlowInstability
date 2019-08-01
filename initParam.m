function paramArray = initParam(P_Pc, Ma, y_pb)
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
    U_ref = Ma * sqrt(P_ref/rho_ref);
    
    height = 0.3e-3; %m
    Rey = U_ref/(mu_ref/height/rho_ref);

    Pe = U_ref*height/d_ref;

    paramArray = {T_base P_ref rho_ref Rgas cp_ref mu_ref lambda_ref Ma Pr Rey h_ref T_ref d_ref Pe};
end

