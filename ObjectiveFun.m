function obj = ObjectiveFun(T, P_ref, x_CH4_ref)
    rho_hat = PC_SAFT_PT_Cubic(P_ref, T, x_CH4_ref, 'CH4', 'O2', 0.5);
    cp_R = PC_SAFT_DT_cp(rho_hat, T, x_CH4_ref, 'CH4', 'O2');
    obj = -cp_R;
end

