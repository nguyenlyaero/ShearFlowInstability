function baseFlowArray = Calc_BaseFlow_Nodes(N, paramArray)
    [xi_vec, ~, ~, ~, ~] = Dmat(N);
    
    [u0__x_vec, du0x_dxi_vec, d2u0x_dxidxi_vec,...
        rho0_vec, drho0_dxi_vec, d2rho0_dxidxi_vec,...
        T0_vec, dT0_dxi_vec, d2T0_dxidxi_vec, ...
        Y0_vec, dY0_dxi_vec, d2Y0_dxidxi_vec, ...
        D12_vec, dD_dT_vec, dD_drho_vec, dD_dY_vec, d2D_dTdT_vec, d2D_dTdrho_vec, d2D_dTdY_vec, d2D_drhodrho_vec, d2D_drhodY_vec, d2D_dYdY_vec,...
        mu_vec, dmu_dT_vec, dmu_drho_vec, dmu_dY_vec, d2mu_dTdT_vec, d2mu_dTdrho_vec, d2mu_dTdY_vec, d2mu_drhodrho_vec, d2mu_drhodY_vec, d2mu_dYdY_vec, ...
        P_vec, dp_dT_vec, dp_drho_vec, dp_dY_vec, d2p_dTdT_vec, d2p_dTdrho_vec, d2p_dTdY_vec, d2p_drhodrho_vec, d2p_drhodY_vec, d2p_dYdY_vec, ...
        h_vec, dh_dT_vec, dh_drho_vec, dh_dY_vec,...
        lambda_vec, dlambda_dT_vec, dlambda_drho_vec, dlambda_dY_vec, d2lambda_dTdT_vec, d2lambda_dTdrho_vec, d2lambda_dTdY_vec, d2lambda_drhodrho_vec, d2lambda_drhodY_vec, d2lambda_dYdY_vec] = deal(zeros(N,1));
    parfor i=1:N
        fprintf('%d \n', i);
         [u0__x_vec(i), du0x_dxi_vec(i), d2u0x_dxidxi_vec(i),...
        rho0_vec(i), drho0_dxi_vec(i), d2rho0_dxidxi_vec(i),...
        T0_vec(i), dT0_dxi_vec(i), d2T0_dxidxi_vec(i), ...
        Y0_vec(i), dY0_dxi_vec(i), d2Y0_dxidxi_vec(i)] = TWO_D_Calc_BaseFlow(xi_vec(i), paramArray);

        [D12_vec(i), dD_dT_vec(i), dD_drho_vec(i), dD_dY_vec(i), d2D_dTdT_vec(i), d2D_dTdrho_vec(i), d2D_dTdY_vec(i), d2D_drhodrho_vec(i), d2D_drhodY_vec(i), d2D_dYdY_vec(i),...
        mu_vec(i), dmu_dT_vec(i), dmu_drho_vec(i), dmu_dY_vec(i), d2mu_dTdT_vec(i), d2mu_dTdrho_vec(i), d2mu_dTdY_vec(i), d2mu_drhodrho_vec(i), d2mu_drhodY_vec(i), d2mu_dYdY_vec(i), ...
        P_vec(i), dp_dT_vec(i), dp_drho_vec(i), dp_dY_vec(i), d2p_dTdT_vec(i), d2p_dTdrho_vec(i), d2p_dTdY_vec(i), d2p_drhodrho_vec(i), d2p_drhodY_vec(i), d2p_dYdY_vec(i), ...
        h_vec(i), dh_dT_vec(i), dh_drho_vec(i), dh_dY_vec(i),...
        lambda_vec(i), dlambda_dT_vec(i), dlambda_drho_vec(i), dlambda_dY_vec(i), d2lambda_dTdT_vec(i), d2lambda_dTdrho_vec(i), d2lambda_dTdY_vec(i), d2lambda_drhodrho_vec(i), d2lambda_drhodY_vec(i), d2lambda_dYdY_vec(i)]...
        = TWO_D_Calc_BaseFlowThermal(T0_vec(i), rho0_vec(i), Y0_vec(i), paramArray);
    end
    baseFlowArray = {u0__x_vec, du0x_dxi_vec, d2u0x_dxidxi_vec,...
        rho0_vec, drho0_dxi_vec, d2rho0_dxidxi_vec,...
        T0_vec, dT0_dxi_vec, d2T0_dxidxi_vec, ...
        Y0_vec, dY0_dxi_vec, d2Y0_dxidxi_vec, ...
        D12_vec, dD_dT_vec, dD_drho_vec, dD_dY_vec, d2D_dTdT_vec, d2D_dTdrho_vec, d2D_dTdY_vec, d2D_drhodrho_vec, d2D_drhodY_vec, d2D_dYdY_vec,...
        mu_vec, dmu_dT_vec, dmu_drho_vec, dmu_dY_vec, d2mu_dTdT_vec, d2mu_dTdrho_vec, d2mu_dTdY_vec, d2mu_drhodrho_vec, d2mu_drhodY_vec, d2mu_dYdY_vec, ...
        P_vec, dp_dT_vec, dp_drho_vec, dp_dY_vec, d2p_dTdT_vec, d2p_dTdrho_vec, d2p_dTdY_vec, d2p_drhodrho_vec, d2p_drhodY_vec, d2p_dYdY_vec, ...
        h_vec, dh_dT_vec, dh_drho_vec, dh_dY_vec,...
        lambda_vec, dlambda_dT_vec, dlambda_drho_vec, dlambda_dY_vec, d2lambda_dTdT_vec, d2lambda_dTdrho_vec, d2lambda_dTdY_vec, d2lambda_drhodrho_vec, d2lambda_drhodY_vec, d2lambda_dYdY_vec};
end

