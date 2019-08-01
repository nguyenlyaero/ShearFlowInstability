function [A,B] = TWO_D_CalcABTemporalCheb(alpha, N, paramArray, baseFlowArray)


[T_base, P_ref, rho_ref, Rgas, cp_ref, mu_ref, lambda_ref, Ma, Pr, Rey, h_ref, T_ref, d_ref, Pe] = deal(paramArray{:});
[u0__x_vec, du0x_dxi_vec, d2u0x_dxidxi_vec,...
        rho0_vec, drho0_dxi_vec, d2rho0_dxidxi_vec,...
        T0_vec, dT0_dxi_vec, d2T0_dxidxi_vec, ...
        Y0_vec, dY0_dxi_vec, d2Y0_dxidxi_vec, ...
        D12_vec, dD_dT_vec, dD_drho_vec, dD_dY_vec, d2D_dTdT_vec, d2D_dTdrho_vec, d2D_dTdY_vec, d2D_drhodrho_vec, d2D_drhodY_vec, d2D_dYdY_vec,...
        mu_vec, dmu_dT_vec, dmu_drho_vec, dmu_dY_vec, d2mu_dTdT_vec, d2mu_dTdrho_vec, d2mu_dTdY_vec, d2mu_drhodrho_vec, d2mu_drhodY_vec, d2mu_dYdY_vec, ...
        P_vec, dp_dT_vec, dp_drho_vec, dp_dY_vec, d2p_dTdT_vec, d2p_dTdrho_vec, d2p_dTdY_vec, d2p_drhodrho_vec, d2p_drhodY_vec, d2p_dYdY_vec, ...
        h_vec, dh_dT_vec, dh_drho_vec, dh_dY_vec,...
        lambda_vec, dlambda_dT_vec, dlambda_drho_vec, dlambda_dY_vec, d2lambda_dTdT_vec, d2lambda_dTdrho_vec, d2lambda_dTdY_vec, d2lambda_drhodrho_vec, d2lambda_drhodY_vec, d2lambda_dYdY_vec] = deal(baseFlowArray{:});

% Load Base Flow and Chebychev polynomials and collocation points

load(['Chebyshev_' sprintf('%d', N) '.mat']);

A = zeros(5*N, 5*N);
B = zeros(5*N, 5*N);
% Equations at interior points
for i=6:5:5*(N-1)
    xi_n = round(i/5)+1;
%     fprintf('%d \n', xi_n);
    xi = xi_vec(xi_n);
    % Functions of xi
    [u0__x, du0x_dxi, d2u0x_dxidxi,...
        rho0, drho0_dxi, d2rho0_dxidxi,...
        T0, dT0_dxi, d2T0_dxidxi, Y0, dY0_dxi, d2Y0_dxidxi] =  deal(u0__x_vec(xi_n), du0x_dxi_vec(xi_n), d2u0x_dxidxi_vec(xi_n),...
        rho0_vec(xi_n), drho0_dxi_vec(xi_n), d2rho0_dxidxi_vec(xi_n),...
        T0_vec(xi_n), dT0_dxi_vec(xi_n), d2T0_dxidxi_vec(xi_n),...
        Y0_vec(xi_n), dY0_dxi_vec(xi_n), d2Y0_dxidxi_vec(xi_n));

    [D12, dD_dT, dD_drho, dD_dY, d2D_dTdT, d2D_dTdrho, d2D_dTdY, d2D_drhodrho, d2D_drhodY, d2D_dYdY,...
        mu, dmu_dT, dmu_drho, dmu_dY, d2mu_dTdT, d2mu_dTdrho, d2mu_dTdY, d2mu_drhodrho, d2mu_drhodY, d2mu_dYdY, ...
        P, dp_dT, dp_drho, dp_dY, d2p_dTdT, d2p_dTdrho, d2p_dTdY, d2p_drhodrho, d2p_drhodY, d2p_dYdY, ...
        h, dh_dT, dh_drho, dh_dY,...
        lambda, dlambda_dT, dlambda_drho, dlambda_dY, d2lambda_dTdT, d2lambda_dTdrho, d2lambda_dTdY, d2lambda_drhodrho, d2lambda_drhodY, d2lambda_dYdY]...
        = deal(D12_vec(xi_n), dD_dT_vec(xi_n), dD_drho_vec(xi_n), dD_dY_vec(xi_n), d2D_dTdT_vec(xi_n), d2D_dTdrho_vec(xi_n), d2D_dTdY_vec(xi_n), d2D_drhodrho_vec(xi_n), d2D_drhodY_vec(xi_n), d2D_dYdY_vec(xi_n),...
        mu_vec(xi_n), dmu_dT_vec(xi_n), dmu_drho_vec(xi_n), dmu_dY_vec(xi_n), d2mu_dTdT_vec(xi_n), d2mu_dTdrho_vec(xi_n), d2mu_dTdY_vec(xi_n), d2mu_drhodrho_vec(xi_n), d2mu_drhodY_vec(xi_n), d2mu_dYdY_vec(xi_n), ...
        P_vec(xi_n), dp_dT_vec(xi_n), dp_drho_vec(xi_n), dp_dY_vec(xi_n), d2p_dTdT_vec(xi_n), d2p_dTdrho_vec(xi_n), d2p_dTdY_vec(xi_n), d2p_drhodrho_vec(xi_n), d2p_drhodY_vec(xi_n), d2p_dYdY_vec(xi_n), ...
        h_vec(xi_n), dh_dT_vec(xi_n), dh_drho_vec(xi_n), dh_dY_vec(xi_n),...
        lambda_vec(xi_n), dlambda_dT_vec(xi_n), dlambda_drho_vec(xi_n), dlambda_dY_vec(xi_n), d2lambda_dTdT_vec(xi_n), d2lambda_dTdrho_vec(xi_n), d2lambda_dTdY_vec(xi_n), d2lambda_drhodrho_vec(xi_n), d2lambda_drhodY_vec(xi_n), d2lambda_dYdY_vec(xi_n));
    % Coefficients for each order of Chebychev polynomial
    for j=1:N
        TCn = D0(xi_n,j);
        dTCn_dxi = D1(xi_n,j);
        d2TCn_dxi2 = D2(xi_n,j);
        
        cg = 1i * TCn * alpha * rho0;
        cg0 = -sqrt((-xi ^ 2 + 1)) * (xi - 1) * (xi + 1) * (TCn * drho0_dxi + rho0 * dTCn_dxi);
        cg1 = 0;
        cg2 = 1i * TCn * alpha * u0__x;
        cg3 = 0;
        cg4 = 0;
        cg5 = rho0 * TCn * dY0_dxi * (-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1);
        cg6 = (xi + 1) ^ 2 * (xi - 1) ^ 2 * (TCn * rho0 * d2D_dTdY * (xi - 1) * (xi + 1) * dY0_dxi ^ 2 + ((((dT0_dxi * d2D_dTdT + drho0_dxi * d2D_dTdrho) * xi ^ 2 + 3 * xi * dD_dT - drho0_dxi * d2D_dTdrho - dT0_dxi * d2D_dTdT) * TCn + dTCn_dxi * dD_dT * (xi - 1) * (xi + 1)) * rho0 + drho0_dxi * TCn * dD_dT * (xi - 1) * (xi + 1)) * dY0_dxi + d2Y0_dxidxi * TCn * rho0 * dD_dT * (xi - 1) * (xi + 1)) / Pe;
        cg7 = (TCn * (xi - 1) * (xi + 1) * (rho0 * d2D_drhodY + dD_dY) * dY0_dxi ^ 2 + ((((dT0_dxi * d2D_dTdrho + d2D_drhodrho * drho0_dxi) * rho0 + 2 * drho0_dxi * dD_drho + dD_dT * dT0_dxi) * xi ^ 2 + (3 * dD_drho * rho0 + 3 * D12) * xi + (-dT0_dxi * d2D_dTdrho - d2D_drhodrho * drho0_dxi) * rho0 - 2 * drho0_dxi * dD_drho - dD_dT * dT0_dxi) * TCn + dTCn_dxi * (xi - 1) * (xi + 1) * (dD_drho * rho0 + D12)) * dY0_dxi + d2Y0_dxidxi * TCn * (xi - 1) * (xi + 1) * (dD_drho * rho0 + D12)) * (xi + 1) ^ 2 * (xi - 1) ^ 2 / Pe;
        cg8 = (((((d2D_dYdY * dY0_dxi ^ 2 + (dT0_dxi * d2D_dTdY + d2D_drhodY * drho0_dxi) * dY0_dxi + d2Y0_dxidxi * dD_dY) * TCn + 2 * dTCn_dxi * dD_dY * dY0_dxi + (dD_dT * dT0_dxi + drho0_dxi * dD_drho) * dTCn_dxi + D12 * d2TCn_dxi2) * xi ^ 6) + ((3 * TCn * dD_dY * dY0_dxi + 3 * D12 * dTCn_dxi) * xi ^ 5) + (((-3 * d2D_dYdY * dY0_dxi ^ 2 + (-3 * dT0_dxi * d2D_dTdY - 3 * d2D_drhodY * drho0_dxi) * dY0_dxi - 3 * d2Y0_dxidxi * dD_dY) * TCn - 6 * dTCn_dxi * dD_dY * dY0_dxi + (-3 * dD_dT * dT0_dxi - 3 * drho0_dxi * dD_drho) * dTCn_dxi - 3 * D12 * d2TCn_dxi2) * xi ^ 4) + ((-6 * TCn * dD_dY * dY0_dxi - 6 * D12 * dTCn_dxi) * xi ^ 3) + (((3 * d2D_dYdY * dY0_dxi ^ 2 + (3 * dT0_dxi * d2D_dTdY + 3 * d2D_drhodY * drho0_dxi) * dY0_dxi + 3 * d2Y0_dxidxi * dD_dY) * TCn + 6 * dTCn_dxi * dD_dY * dY0_dxi + (3 * dD_dT * dT0_dxi + 3 * drho0_dxi * dD_drho) * dTCn_dxi + 3 * D12 * d2TCn_dxi2) * xi ^ 2) + ((3 * TCn * dD_dY * dY0_dxi + 3 * D12 * dTCn_dxi) * xi) + (-(d2D_dYdY * dY0_dxi ^ 2) + ((-dT0_dxi * d2D_dTdY - d2D_drhodY * drho0_dxi) * dY0_dxi) + 1i * u0__x * alpha * Pe + D12 * alpha ^ 2 - (d2Y0_dxidxi * dD_dY)) * TCn - (2 * dTCn_dxi * dD_dY * dY0_dxi) + ((-dD_dT * dT0_dxi - drho0_dxi * dD_drho) * dTCn_dxi) - (D12 * d2TCn_dxi2)) * rho0 + (drho0_dxi * (xi - 1) ^ 3 * (xi + 1) ^ 3 * (TCn * dD_dY * dY0_dxi + D12 * dTCn_dxi))) / Pe;
        cg9 = ((0.9e1 * ((xi + 1) ^ 2) * ((xi - 1) ^ 2) * ((dmu_dT * dT0_dxi / 0.3e1 + dmu_dY * dY0_dxi / 0.3e1 + dmu_drho * drho0_dxi / 0.3e1) * (xi ^ 2) + (mu * xi) - dmu_dT * dT0_dxi / 0.3e1 - dmu_dY * dY0_dxi / 0.3e1 - dmu_drho * drho0_dxi / 0.3e1) * dTCn_dxi) + (3 * mu * xi ^ 6 * d2TCn_dxi2) - (9 * mu * xi ^ 4 * d2TCn_dxi2) + (9 * mu * xi ^ 2 * d2TCn_dxi2) + ((4 * TCn * alpha ^ 2 - 3 * d2TCn_dxi2) * mu) + 3*1i * rho0 * u0__x * TCn * alpha * Rey) / Rey / 3;
        cg10 = -(-(3 * rho0 * du0x_dxi * TCn * Rey) + 1i * alpha * ((3 * dmu_dT * dT0_dxi + 3 * dmu_dY * dY0_dxi + 3 * dmu_drho * drho0_dxi) * TCn + mu * dTCn_dxi)) * ((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) / Rey / 3;
        cg11 = (((xi + 1) ^ 2 * ((((d2mu_dTdT * dT0_dxi + d2mu_dTdY * dY0_dxi + drho0_dxi * d2mu_dTdrho) * xi ^ 2 + 3 * dmu_dT * xi - d2mu_dTdT * dT0_dxi - d2mu_dTdY * dY0_dxi - drho0_dxi * d2mu_dTdrho) * du0x_dxi + d2u0x_dxidxi * dmu_dT * (xi - 1) * (xi + 1)) * TCn + dTCn_dxi * du0x_dxi * dmu_dT * (xi - 1) * (xi + 1)) * (xi - 1) ^ 2 * Ma ^ 2) + 1i * dp_dT * TCn * alpha * Rey) / (Ma ^ 2) / Rey;
        cg12 = (((xi + 1) ^ 2 * (xi - 1) ^ 2 * ((((d2mu_dTdrho * dT0_dxi + d2mu_drhodY * dY0_dxi + drho0_dxi * d2mu_drhodrho) * xi ^ 2 + 3 * dmu_drho * xi - drho0_dxi * d2mu_drhodrho - d2mu_dTdrho * dT0_dxi - d2mu_drhodY * dY0_dxi) * du0x_dxi + dmu_drho * d2u0x_dxidxi * (xi - 1) * (xi + 1)) * TCn + du0x_dxi * dmu_drho * dTCn_dxi * (xi - 1) * (xi + 1)) * Ma ^ 2) + 1i * dp_drho * TCn * alpha * Rey) / (Ma ^ 2) / Rey;
        cg13 = (((xi + 1) ^ 2 * ((((d2mu_dTdY * dT0_dxi + d2mu_dYdY * dY0_dxi + drho0_dxi * d2mu_drhodY) * xi ^ 2 + 3 * dmu_dY * xi - d2mu_dTdY * dT0_dxi - d2mu_dYdY * dY0_dxi - drho0_dxi * d2mu_drhodY) * du0x_dxi + d2u0x_dxidxi * dmu_dY * (xi - 1) * (xi + 1)) * TCn + dTCn_dxi * du0x_dxi * dmu_dY * (xi - 1) * (xi + 1)) * (xi - 1) ^ 2 * Ma ^ 2) + 1i * dp_dY * TCn * alpha * Rey) / (Ma ^ 2) / Rey;
        cg14 = -0.1e1 / 0.3e1*1i * alpha * ((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) * ((-2 * dmu_dT * dT0_dxi - 2 * dmu_dY * dY0_dxi - 2 * dmu_drho * drho0_dxi) * TCn + mu * dTCn_dxi) / Rey;
        cg15 = ((0.12e2 * ((xi + 1) ^ 2) * ((xi - 1) ^ 2) * ((dmu_dT * dT0_dxi / 0.3e1 + dmu_dY * dY0_dxi / 0.3e1 + dmu_drho * drho0_dxi / 0.3e1) * (xi ^ 2) + (mu * xi) - dmu_dT * dT0_dxi / 0.3e1 - dmu_dY * dY0_dxi / 0.3e1 - dmu_drho * drho0_dxi / 0.3e1) * dTCn_dxi) + (4 * mu * xi ^ 6 * d2TCn_dxi2) - (12 * mu * xi ^ 4 * d2TCn_dxi2) + (12 * mu * xi ^ 2 * d2TCn_dxi2) + ((3 * TCn * alpha ^ 2 - 4 * d2TCn_dxi2) * mu) + 3*1i * rho0 * u0__x * TCn * alpha * Rey) / Rey / 3;
        cg16 = -((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) * (((-dT0_dxi * d2p_dTdT - dY0_dxi * d2p_dTdY - d2p_dTdrho * drho0_dxi) * TCn - dp_dT * dTCn_dxi) * Rey + 1i * Ma ^ 2 * TCn * alpha * dmu_dT * du0x_dxi) / Ma ^ 2 / Rey;
        cg17 = -((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) * (((-dT0_dxi * d2p_dTdrho - dY0_dxi * d2p_drhodY - d2p_drhodrho * drho0_dxi) * TCn - dp_drho * dTCn_dxi) * Rey + 1i * Ma ^ 2 * TCn * alpha * dmu_drho * du0x_dxi) / Ma ^ 2 / Rey;
        cg18 = -((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) * (((-dT0_dxi * d2p_dTdY - dY0_dxi * d2p_dYdY - d2p_drhodY * drho0_dxi) * TCn - dp_dY * dTCn_dxi) * Rey + 1i * Ma ^ 2 * TCn * alpha * dmu_dY * du0x_dxi) / Ma ^ 2 / Rey;
        cg19 = Ma ^ 2 * (((((dmu_dT * dT0_dxi + dmu_dY * dY0_dxi + dmu_drho * drho0_dxi) * dTCn_dxi + mu * d2TCn_dxi2) * u0__x + ((dmu_dT * dT0_dxi + dmu_dY * dY0_dxi + dmu_drho * drho0_dxi) * du0x_dxi + mu * d2u0x_dxidxi) * TCn + 2 * du0x_dxi * mu * dTCn_dxi) * xi ^ 6) + (3 * mu * (TCn * du0x_dxi + u0__x * dTCn_dxi) * xi ^ 5) + ((((-3 * dmu_dT * dT0_dxi - 3 * dmu_dY * dY0_dxi - 3 * dmu_drho * drho0_dxi) * dTCn_dxi - 3 * mu * d2TCn_dxi2) * u0__x + ((-3 * dmu_dT * dT0_dxi - 3 * dmu_dY * dY0_dxi - 3 * dmu_drho * drho0_dxi) * du0x_dxi - 3 * mu * d2u0x_dxidxi) * TCn - 6 * du0x_dxi * mu * dTCn_dxi) * xi ^ 4) - (6 * mu * (TCn * du0x_dxi + u0__x * dTCn_dxi) * xi ^ 3) + ((((3 * dmu_dT * dT0_dxi + 3 * dmu_dY * dY0_dxi + 3 * dmu_drho * drho0_dxi) * dTCn_dxi + 3 * mu * d2TCn_dxi2) * u0__x + ((3 * dmu_dT * dT0_dxi + 3 * dmu_dY * dY0_dxi + 3 * dmu_drho * drho0_dxi) * du0x_dxi + 3 * mu * d2u0x_dxidxi) * TCn + 6 * du0x_dxi * mu * dTCn_dxi) * xi ^ 2) + (3 * mu * (TCn * du0x_dxi + u0__x * dTCn_dxi) * xi) + 1i * rho0 * TCn * alpha * (u0__x ^ 2) * Rey + ((0.4e1 / 0.3e1) * mu * TCn * alpha ^ 2 + ((-dmu_dT * dT0_dxi - dmu_dY * dY0_dxi - dmu_drho * drho0_dxi) * dTCn_dxi) - (mu * d2TCn_dxi2)) * u0__x + (((-dmu_dT * dT0_dxi - dmu_dY * dY0_dxi - dmu_drho * drho0_dxi) * du0x_dxi - mu * d2u0x_dxidxi) * TCn) - (2 * du0x_dxi * mu * dTCn_dxi)) * Rgas / cp_ref / Rey;
        cg20 = -((-xi ^ 2 + 1) ^ (0.3e1 / 0.2e1)) * (Ma ^ 2 * Rgas * (-(rho0 * du0x_dxi * TCn * Rey) + 1i * alpha * (dmu_dT * dT0_dxi * TCn + dmu_dY * dY0_dxi * TCn + dmu_drho * drho0_dxi * TCn + mu * dTCn_dxi / 0.3e1)) * u0__x + 2 * TCn * (-(cp_ref * Rey * rho0 * (dT0_dxi * dh_dT + dY0_dxi * dh_dY + dh_drho * drho0_dxi) / 0.2e1) + 1i * Rgas * Ma ^ 2 * alpha * mu * du0x_dxi)) / cp_ref / Rey;
        cg21 = (((((dT0_dxi ^ 2 * d2lambda_dTdT + (dY0_dxi * d2lambda_dTdY + d2lambda_dTdrho * drho0_dxi) * dT0_dxi + dlambda_dT * d2T0_dxidxi) * cp_ref + Ma ^ 2 * Pr * Rgas * (((d2mu_dTdT * dT0_dxi + d2mu_dTdY * dY0_dxi + drho0_dxi * d2mu_dTdrho) * du0x_dxi + d2u0x_dxidxi * dmu_dT) * u0__x + dmu_dT * du0x_dxi ^ 2)) * TCn + (2 * dT0_dxi * dlambda_dT * dTCn_dxi + (dY0_dxi * dlambda_dY + dlambda_drho * drho0_dxi) * dTCn_dxi + d2TCn_dxi2 * lambda) * cp_ref + Rgas * Ma ^ 2 * du0x_dxi * u0__x * dmu_dT * dTCn_dxi * Pr) * xi ^ 6) + (((3 * Ma ^ 2 * Pr * Rgas * dmu_dT * u0__x * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_dT) * TCn + 3 * cp_ref * lambda * dTCn_dxi) * xi ^ 5) + ((((-3 * dT0_dxi ^ 2 * d2lambda_dTdT + (-3 * dY0_dxi * d2lambda_dTdY - 3 * d2lambda_dTdrho * drho0_dxi) * dT0_dxi - 3 * dlambda_dT * d2T0_dxidxi) * cp_ref - 3 * Ma ^ 2 * Pr * Rgas * (((d2mu_dTdT * dT0_dxi + d2mu_dTdY * dY0_dxi + drho0_dxi * d2mu_dTdrho) * du0x_dxi + d2u0x_dxidxi * dmu_dT) * u0__x + dmu_dT * du0x_dxi ^ 2)) * TCn + (-6 * dT0_dxi * dlambda_dT * dTCn_dxi + (-3 * dY0_dxi * dlambda_dY - 3 * dlambda_drho * drho0_dxi) * dTCn_dxi - 3 * d2TCn_dxi2 * lambda) * cp_ref - 3 * Rgas * Ma ^ 2 * du0x_dxi * u0__x * dmu_dT * dTCn_dxi * Pr) * xi ^ 4) + (((-6 * Ma ^ 2 * Pr * Rgas * dmu_dT * u0__x * du0x_dxi - 6 * cp_ref * dT0_dxi * dlambda_dT) * TCn - 6 * cp_ref * lambda * dTCn_dxi) * xi ^ 3) + ((((3 * dT0_dxi ^ 2 * d2lambda_dTdT + (3 * dY0_dxi * d2lambda_dTdY + 3 * d2lambda_dTdrho * drho0_dxi) * dT0_dxi + 3 * dlambda_dT * d2T0_dxidxi) * cp_ref + 3 * Ma ^ 2 * Pr * Rgas * (((d2mu_dTdT * dT0_dxi + d2mu_dTdY * dY0_dxi + drho0_dxi * d2mu_dTdrho) * du0x_dxi + d2u0x_dxidxi * dmu_dT) * u0__x + dmu_dT * du0x_dxi ^ 2)) * TCn + (6 * dT0_dxi * dlambda_dT * dTCn_dxi + (3 * dY0_dxi * dlambda_dY + 3 * dlambda_drho * drho0_dxi) * dTCn_dxi + 3 * d2TCn_dxi2 * lambda) * cp_ref + 3 * Rgas * Ma ^ 2 * du0x_dxi * u0__x * dmu_dT * dTCn_dxi * Pr) * xi ^ 2) + (((3 * Ma ^ 2 * Pr * Rgas * dmu_dT * u0__x * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_dT) * TCn + 3 * cp_ref * lambda * dTCn_dxi) * xi) + ((1i * rho0 * u0__x * dh_dT * alpha * Pr * Rey - (dT0_dxi ^ 2 * d2lambda_dTdT) + ((-dY0_dxi * d2lambda_dTdY - d2lambda_dTdrho * drho0_dxi) * dT0_dxi) + lambda * alpha ^ 2 - (dlambda_dT * d2T0_dxidxi)) * cp_ref - (Ma ^ 2 * Pr * Rgas * (((d2mu_dTdT * dT0_dxi + d2mu_dTdY * dY0_dxi + drho0_dxi * d2mu_dTdrho) * du0x_dxi + d2u0x_dxidxi * dmu_dT) * u0__x + dmu_dT * du0x_dxi ^ 2))) * TCn + ((-2 * dT0_dxi * dlambda_dT * dTCn_dxi + (-dY0_dxi * dlambda_dY - dlambda_drho * drho0_dxi) * dTCn_dxi - d2TCn_dxi2 * lambda) * cp_ref) - (Rgas * Ma ^ 2 * du0x_dxi * u0__x * dmu_dT * dTCn_dxi * Pr)) / Pr / Rey / cp_ref;
        cg22 = ((((Ma ^ 2 * (((d2mu_dTdrho * dT0_dxi + d2mu_drhodY * dY0_dxi + drho0_dxi * d2mu_drhodrho) * du0x_dxi + dmu_drho * d2u0x_dxidxi) * u0__x + dmu_drho * du0x_dxi ^ 2) * Rgas * Pr + cp_ref * (dT0_dxi ^ 2 * d2lambda_dTdrho + (dY0_dxi * d2lambda_drhodY + drho0_dxi * d2lambda_drhodrho) * dT0_dxi + dlambda_drho * d2T0_dxidxi)) * xi ^ 6) + ((3 * Ma ^ 2 * Pr * Rgas * u0__x * dmu_drho * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_drho) * xi ^ 5) + ((-3 * Ma ^ 2 * (((d2mu_dTdrho * dT0_dxi + d2mu_drhodY * dY0_dxi + drho0_dxi * d2mu_drhodrho) * du0x_dxi + dmu_drho * d2u0x_dxidxi) * u0__x + dmu_drho * du0x_dxi ^ 2) * Rgas * Pr - 3 * cp_ref * (dT0_dxi ^ 2 * d2lambda_dTdrho + (dY0_dxi * d2lambda_drhodY + drho0_dxi * d2lambda_drhodrho) * dT0_dxi + dlambda_drho * d2T0_dxidxi)) * xi ^ 4) + ((-6 * Ma ^ 2 * Pr * Rgas * u0__x * dmu_drho * du0x_dxi - 6 * cp_ref * dT0_dxi * dlambda_drho) * xi ^ 3) + ((3 * Ma ^ 2 * (((d2mu_dTdrho * dT0_dxi + d2mu_drhodY * dY0_dxi + drho0_dxi * d2mu_drhodrho) * du0x_dxi + dmu_drho * d2u0x_dxidxi) * u0__x + dmu_drho * du0x_dxi ^ 2) * Rgas * Pr + 3 * cp_ref * (dT0_dxi ^ 2 * d2lambda_dTdrho + (dY0_dxi * d2lambda_drhodY + drho0_dxi * d2lambda_drhodrho) * dT0_dxi + dlambda_drho * d2T0_dxidxi)) * xi ^ 2) + ((3 * Ma ^ 2 * Pr * Rgas * u0__x * dmu_drho * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_drho) * xi) + (-(Ma ^ 2 * (((d2mu_dTdrho * dT0_dxi + d2mu_drhodY * dY0_dxi + drho0_dxi * d2mu_drhodrho) * du0x_dxi + dmu_drho * d2u0x_dxidxi) * u0__x + dmu_drho * du0x_dxi ^ 2) * Rgas) + 1i * rho0 * u0__x * dh_drho * alpha * Rey * cp_ref) * Pr - (cp_ref * (dT0_dxi ^ 2 * d2lambda_dTdrho + (dY0_dxi * d2lambda_drhodY + drho0_dxi * d2lambda_drhodrho) * dT0_dxi + dlambda_drho * d2T0_dxidxi))) * TCn + (dTCn_dxi * (xi - 1) ^ 3 * (xi + 1) ^ 3 * (Ma ^ 2 * Pr * Rgas * u0__x * dmu_drho * du0x_dxi + cp_ref * dT0_dxi * dlambda_drho))) / Pr / Rey / cp_ref;
        cg23 = ((((Ma ^ 2 * (((d2mu_dTdY * dT0_dxi + d2mu_dYdY * dY0_dxi + drho0_dxi * d2mu_drhodY) * du0x_dxi + d2u0x_dxidxi * dmu_dY) * u0__x + dmu_dY * du0x_dxi ^ 2) * Rgas * Pr + (dT0_dxi ^ 2 * d2lambda_dTdY + (dY0_dxi * d2lambda_dYdY + d2lambda_drhodY * drho0_dxi) * dT0_dxi + dlambda_dY * d2T0_dxidxi) * cp_ref) * xi ^ 6) + ((3 * Ma ^ 2 * Pr * Rgas * dmu_dY * u0__x * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_dY) * xi ^ 5) + ((-3 * Ma ^ 2 * (((d2mu_dTdY * dT0_dxi + d2mu_dYdY * dY0_dxi + drho0_dxi * d2mu_drhodY) * du0x_dxi + d2u0x_dxidxi * dmu_dY) * u0__x + dmu_dY * du0x_dxi ^ 2) * Rgas * Pr - 3 * (dT0_dxi ^ 2 * d2lambda_dTdY + (dY0_dxi * d2lambda_dYdY + d2lambda_drhodY * drho0_dxi) * dT0_dxi + dlambda_dY * d2T0_dxidxi) * cp_ref) * xi ^ 4) + ((-6 * Ma ^ 2 * Pr * Rgas * dmu_dY * u0__x * du0x_dxi - 6 * cp_ref * dT0_dxi * dlambda_dY) * xi ^ 3) + ((3 * Ma ^ 2 * (((d2mu_dTdY * dT0_dxi + d2mu_dYdY * dY0_dxi + drho0_dxi * d2mu_drhodY) * du0x_dxi + d2u0x_dxidxi * dmu_dY) * u0__x + dmu_dY * du0x_dxi ^ 2) * Rgas * Pr + 3 * (dT0_dxi ^ 2 * d2lambda_dTdY + (dY0_dxi * d2lambda_dYdY + d2lambda_drhodY * drho0_dxi) * dT0_dxi + dlambda_dY * d2T0_dxidxi) * cp_ref) * xi ^ 2) + ((3 * Ma ^ 2 * Pr * Rgas * dmu_dY * u0__x * du0x_dxi + 3 * cp_ref * dT0_dxi * dlambda_dY) * xi) + (-(Ma ^ 2 * (((d2mu_dTdY * dT0_dxi + d2mu_dYdY * dY0_dxi + drho0_dxi * d2mu_drhodY) * du0x_dxi + d2u0x_dxidxi * dmu_dY) * u0__x + dmu_dY * du0x_dxi ^ 2) * Rgas) + 1i * rho0 * u0__x * dh_dY * alpha * Rey * cp_ref) * Pr - ((dT0_dxi ^ 2 * d2lambda_dTdY + (dY0_dxi * d2lambda_dYdY + d2lambda_drhodY * drho0_dxi) * dT0_dxi + dlambda_dY * d2T0_dxidxi) * cp_ref)) * TCn + (dTCn_dxi * (xi - 1) ^ 3 * (xi + 1) ^ 3 * (Ma ^ 2 * Pr * Rgas * dmu_dY * u0__x * du0x_dxi + cp_ref * dT0_dxi * dlambda_dY))) / Pr / Rey / cp_ref;




        A(i,j) = cg;
        A(i,N+j) = cg0;
        A(i,2*N+j) = cg1;
        A(i,3*N+j) = cg2;
        A(i,4*N+j) = cg3;
        
        A(i+1,j) = cg4;
        A(i+1,N+j) = cg5;
        A(i+1,2*N+j) = cg6;
        A(i+1,3*N+j) = cg7;
        A(i+1,4*N+j) = cg8;

        A(i+2,j) = cg9;
        A(i+2,N+j) = cg10;
        A(i+2,2*N+j) = cg11;
        A(i+2,3*N+j) = cg12;
        A(i+2,4*N+j) = cg13;

        A(i+3,j) = cg14;
        A(i+3,N+j) = cg15;
        A(i+3,2*N+j) = cg16;
        A(i+3,3*N+j) = cg17;
        A(i+3,4*N+j) = cg18;
        
        A(i+4,j) = cg19;
        A(i+4,N+j) = cg20;
        A(i+4,2*N+j) = cg21;
        A(i+4,3*N+j) = cg22;
        A(i+4,4*N+j) = cg23;
        
        cg24 = 0;
        cg25 = 0;
        cg26 = 0;
        cg27 = 1i * TCn;
        cg28 = 0;
        cg29 = 0;
        cg30 = 0;
        cg31 = 0;
        cg32 = 0;
        cg33 = 1i * rho0 * TCn;
        cg34 = 1i * rho0 * TCn;
        cg35 = 0;
        cg36 = 0;
        cg37 = 0;
        cg38 = 0;
        cg39 = 0;
        cg40 = 1i * rho0 * TCn;
        cg41 = 0;
        cg42 = 0;
        cg43 = 0;
        cg44 = 1i * rho0 * Rgas / cp_ref * Ma ^ 2 * u0__x * TCn;
        cg45 = 0;
        cg46 = -1i * TCn * (-cp_ref * dh_dT * rho0 + Rgas * dp_dT) / cp_ref;
        cg47 = -1i * TCn * (-cp_ref * dh_drho * rho0 + Rgas * dp_drho) / cp_ref;
        cg48 = -1i * TCn * (-cp_ref * dh_dY * rho0 + Rgas * dp_dY) / cp_ref;
        
        B(i,j) = cg24;
        B(i,N+j) = cg25;
        B(i,2*N+j) = cg26;
        B(i,3*N+j) = cg27;
        B(i,4*N+j) = cg28;
        
        B(i+1,j) = cg29;
        B(i+1,N+j) = cg30;
        B(i+1,2*N+j) = cg31;
        B(i+1,3*N+j) = cg32;
        B(i+1,4*N+j) = cg33;

        B(i+2,j) = cg34;
        B(i+2,N+j) = cg35;
        B(i+2,2*N+j) = cg36;
        B(i+2,3*N+j) = cg37;
        B(i+2,4*N+j) = cg38;

        B(i+3,j) = cg39;
        B(i+3,N+j) = cg40;
        B(i+3,2*N+j) = cg41;
        B(i+3,3*N+j) = cg42;
        B(i+3,4*N+j) = cg43;
        
        B(i+4,j) = cg44;
        B(i+4,N+j) = cg45;
        B(i+4,2*N+j) = cg46;
        B(i+4,3*N+j) = cg47;
        B(i+4,4*N+j) = cg48;
    end  
end

% Boundary Conditions
C=-complex(500,500);

for i=1:5
    for j=1:N
        % Bottom
        A(i,N*(i-1)+j) = C*D0(1,j);
        B(i,N*(i-1)+j) = D0(1,j);
        % Top
        A(5*N-(i-1),N*(i-1)+j) = C*D0(end,j);
        B(5*N-(i-1),N*(i-1)+j) = D0(end,j);
    end
end

end