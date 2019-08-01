function [u0__x, du0x_dxi, d2u0x_dxidxi,...
        rho0, drho0_dxi, d2rho0_dxidxi,...
        T0, dT0_dxi, d2T0_dxidxi, ...
        Y0, dY0_dxi, d2Y0_dxidxi] = TWO_D_Calc_BaseFlow(xi, paramArray)
    [T_base, P_ref, rho_ref, Rgas, cp_ref, mu_ref, lambda_ref, Ma, Pr, Rey, h_ref, T_ref, d_ref, Pe] = deal(paramArray{:});
    
    diff1 = ((-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1) + xi ^ 2 * (-xi ^ 2 + 1) ^ (-0.3e1 / 0.2e1)) * (0.1e1 - tanh((xi * (-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1))) ^ 2);
    diff2 = (3 * (-xi ^ 2 + 1) ^ (-0.3e1 / 0.2e1) * xi + 3 * xi ^ 3 * (-xi ^ 2 + 1) ^ (-0.5e1 / 0.2e1)) * (0.1e1 - tanh((xi * (-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1))) ^ 2) - 0.2e1 * (((-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1) + xi ^ 2 * (-xi ^ 2 + 1) ^ (-0.3e1 / 0.2e1)) ^ 2) * tanh((xi * (-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1))) * (0.1e1 - tanh((xi * (-xi ^ 2 + 1) ^ (-0.1e1 / 0.2e1))) ^ 2);
    
    u0__x = 1+3/4*tanh(xi/sqrt(1-xi^2));
    if xi ~= 1 && xi ~=-1
        du0x_dxi = 3/4*diff1;
        d2u0x_dxidxi = 3/4*diff2;
    else
        du0x_dxi = 0;
        d2u0x_dxidxi = 0;
    end
    
    T0 = T_base/T_ref*(1+0.2*tanh(xi/sqrt(1-xi^2)));
    if xi ~= 1 && xi ~=-1
        dT0_dxi = 0.2*diff1;
        dT0_dxi = T_base/T_ref* dT0_dxi;
        d2T0_dxidxi = 0.2*diff2;
        d2T0_dxidxi = T_base/T_ref*d2T0_dxidxi;
    else
        dT0_dxi = 0;
        d2T0_dxidxi = 0;
    end
    
    Y0 = 0.5*(1+tanh(xi/sqrt(1-xi^2)));
    if xi ~= 1 && xi ~=-1
        dY0_dxi = 0.5*diff1;
        d2Y0_dxidxi = 0.5*diff2;
    else
        dY0_dxi = 0;
        d2Y0_dxidxi = 0;
    end
    
    rho0 = refprop_TP('rho', T0, 1, Y0, paramArray);
    
    grad_rho = gradest(@(vec) refprop_TP('rho', vec(1), 1, vec(2), paramArray), [T0 Y0]);
    hess_rho = hessian(@(vec) refprop_TP('rho', vec(1), 1, vec(2), paramArray), [T0 Y0]);
    
    drho0_dT = grad_rho(1);
    drho0_dY = grad_rho(2);
    d2rho0_dTdT = hess_rho(1,1);
    d2rho0_dYdY = hess_rho(2,2);
    
    drho0_dxi = drho0_dT*dT0_dxi + drho0_dY*dY0_dxi;
    d2rho0_dxidxi = drho0_dT*d2T0_dxidxi + d2rho0_dTdT*(dT0_dxi)^2 + drho0_dY*d2Y0_dxidxi + d2rho0_dYdY*(dY0_dxi)^2;
   
end