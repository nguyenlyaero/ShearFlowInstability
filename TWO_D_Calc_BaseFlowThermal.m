function [D12, dD_dT, dD_drho, dD_dY, d2D_dTdT, d2D_dTdrho, d2D_dTdY, d2D_drhodrho, d2D_drhodY, d2D_dYdY,...
        mu, dmu_dT, dmu_drho, dmu_dY, d2mu_dTdT, d2mu_dTdrho, d2mu_dTdY, d2mu_drhodrho, d2mu_drhodY, d2mu_dYdY, ...
        P, dp_dT, dp_drho, dp_dY, d2p_dTdT, d2p_dTdrho, d2p_dTdY, d2p_drhodrho, d2p_drhodY, d2p_dYdY, ...
        h, dh_dT, dh_drho, dh_dY,...
        lambda, dlambda_dT, dlambda_drho, dlambda_dY, d2lambda_dTdT, d2lambda_dTdrho, d2lambda_dTdY, d2lambda_drhodrho, d2lambda_drhodY, d2lambda_dYdY]...
        = TWO_D_Calc_BaseFlowThermal(T0, rho0, Y0, paramArray)
    
%     step = 1e-6;
%     
%     D12 = refprop_TD('d', T0, rho0, Y0);
%     dD_dT = (refprop_TD('d', T0*(1+step), rho0, Y0) - refprop_TD('d', T0*(1-step), rho0, Y0))/(2*step*T0);
%     dD_drho = (refprop_TD('d', T0, rho0*(1+step), Y0) - refprop_TD('d', T0, rho0*(1-step), Y0))/(2*step*rho0);
%     dD_dY = (refprop_TD('d', T0, rho0, Y0*(1+step)) - refprop_TD('d', T0, rho0, Y0*(1-step)))/(2*step*Y0);
%     d2D_dTdT = (refprop_TD('d', T0*(1+step), rho0, Y0) - 2*D12 + refprop_TD('d', T0*(1-step), rho0, Y0))/(step*T0)^2;
%     d2D_dTdrho = (refprop_TD('d', T0*(1+step), rho0*(1+step), Y0)...
%         - refprop_TD('d', T0*(1+step), rho0*(1-step), Y0)...
%         - refprop_TD('d', T0*(1-step), rho0*(1+step), Y0)...
%         + refprop_TD('d', T0*(1-step), rho0*(1-step), Y0))/(4*step^2*T0*rho0);
%     d2D_dTdY = (refprop_TD('d', T0*(1+step), rho0, Y0*(1+step))...
%         - refprop_TD('d', T0*(1+step), rho0, Y0*(1-step))...
%         - refprop_TD('d', T0*(1-step), rho0, Y0*(1+step))...
%         + refprop_TD('d', T0*(1-step), rho0, Y0*(1-step)))/(4*step^2*T0*Y0);
%     d2D_drhodrho = (refprop_TD('d', T0, rho0*(1+step), Y0) - 2*D12 + refprop_TD('d', T0, rho0*(1-step), Y0))/(step*rho0)^2;
%     d2D_drhodY = (refprop_TD('d', T0, rho0*(1+step), Y0*(1+step))...
%         - refprop_TD('d', T0, rho0*(1+step), Y0*(1-step))...
%         - refprop_TD('d', T0, rho0*(1-step), Y0*(1+step))...
%         + refprop_TD('d', T0, rho0*(1-step), Y0*(1-step)))/(4*step^2*rho0*Y0);
%     d2D_dYdY = (refprop_TD('d', T0, rho0, Y0*(1+step)) - 2*D12 + refprop_TD('d', T0, rho0, Y0*(1-step)))/(step*Y0)^2;
%     
%     
%     mu = refprop_TD('mu', T0, rho0, Y0);
%     dmu_dT = (refprop_TD('mu', T0*(1+step), rho0, Y0) - refprop_TD('mu', T0*(1-step), rho0, Y0))/(2*step*T0);
%     dmu_drho = (refprop_TD('mu', T0, rho0*(1+step), Y0) - refprop_TD('mu', T0, rho0*(1-step), Y0))/(2*step*rho0);
%     dmu_dY = (refprop_TD('mu', T0, rho0, Y0*(1+step)) - refprop_TD('mu', T0, rho0, Y0*(1-step)))/(2*step*Y0);
%     d2mu_dTdT = (refprop_TD('mu', T0*(1+step), rho0, Y0) - 2*D12 + refprop_TD('mu', T0*(1-step), rho0, Y0))/(step*T0)^2;
%     d2mu_dTdrho = (refprop_TD('mu', T0*(1+step), rho0*(1+step), Y0)...
%         - refprop_TD('mu', T0*(1+step), rho0*(1-step), Y0)...
%         - refprop_TD('mu', T0*(1-step), rho0*(1+step), Y0)...
%         + refprop_TD('mu', T0*(1-step), rho0*(1-step), Y0))/(4*step^2*T0*rho0);
%     d2mu_dTdY = (refprop_TD('mu', T0*(1+step), rho0, Y0*(1+step))...
%         - refprop_TD('mu', T0*(1+step), rho0, Y0*(1-step))...
%         - refprop_TD('mu', T0*(1-step), rho0, Y0*(1+step))...
%         + refprop_TD('mu', T0*(1-step), rho0, Y0*(1-step)))/(4*step^2*T0*Y0);
%     d2mu_drhodrho = (refprop_TD('mu', T0, rho0*(1+step), Y0) - 2*D12 + refprop_TD('mu', T0, rho0*(1-step), Y0))/(step*rho0)^2;
%     d2mu_drhodY = (refprop_TD('mu', T0, rho0*(1+step), Y0*(1+step))...
%         - refprop_TD('mu', T0, rho0*(1+step), Y0*(1-step))...
%         - refprop_TD('mu', T0, rho0*(1-step), Y0*(1+step))...
%         + refprop_TD('mu', T0, rho0*(1-step), Y0*(1-step)))/(4*step^2*rho0*Y0);
%     d2mu_dYdY = (refprop_TD('mu', T0, rho0, Y0*(1+step)) - 2*D12 + refprop_TD('mu', T0, rho0, Y0*(1-step)))/(step*Y0)^2;
%     
%     P = refprop_TD('P', T0, rho0, Y0);
%     dp_dT = (refprop_TD('P', T0*(1+step), rho0, Y0) - refprop_TD('P', T0*(1-step), rho0, Y0))/(2*step*T0);
%     dp_drho = (refprop_TD('P', T0, rho0*(1+step), Y0) - refprop_TD('P', T0, rho0*(1-step), Y0))/(2*step*rho0);
%     dp_dY = (refprop_TD('P', T0, rho0, Y0*(1+step)) - refprop_TD('P', T0, rho0, Y0*(1-step)))/(2*step*Y0);
%     d2p_dTdT = (refprop_TD('P', T0*(1+step), rho0, Y0) - 2*D12 + refprop_TD('P', T0*(1-step), rho0, Y0))/(step*T0)^2;
%     d2p_dTdrho = (refprop_TD('P', T0*(1+step), rho0*(1+step), Y0)...
%         - refprop_TD('P', T0*(1+step), rho0*(1-step), Y0)...
%         - refprop_TD('P', T0*(1-step), rho0*(1+step), Y0)...
%         + refprop_TD('P', T0*(1-step), rho0*(1-step), Y0))/(4*step^2*T0*rho0);
%     d2p_dTdY = (refprop_TD('P', T0*(1+step), rho0, Y0*(1+step))...
%         - refprop_TD('P', T0*(1+step), rho0, Y0*(1-step))...
%         - refprop_TD('P', T0*(1-step), rho0, Y0*(1+step))...
%         + refprop_TD('P', T0*(1-step), rho0, Y0*(1-step)))/(4*step^2*T0*Y0);
%     d2p_drhodrho = (refprop_TD('P', T0, rho0*(1+step), Y0) - 2*D12 + refprop_TD('P', T0, rho0*(1-step), Y0))/(step*rho0)^2;
%     d2p_drhodY = (refprop_TD('P', T0, rho0*(1+step), Y0*(1+step))...
%         - refprop_TD('P', T0, rho0*(1+step), Y0*(1-step))...
%         - refprop_TD('P', T0, rho0*(1-step), Y0*(1+step))...
%         + refprop_TD('P', T0, rho0*(1-step), Y0*(1-step)))/(4*step^2*rho0*Y0);
%     d2p_dYdY = (refprop_TD('P', T0, rho0, Y0*(1+step)) - 2*D12 + refprop_TD('P', T0, rho0, Y0*(1-step)))/(step*Y0)^2;
% 
%     h = refprop_TD('h', T0, rho0, Y0);
%     dh_dT = (refprop_TD('h', T0*(1+step), rho0, Y0) - refprop_TD('h', T0*(1-step), rho0, Y0))/(2*step*T0);
%     dh_drho = (refprop_TD('h', T0, rho0*(1+step), Y0) - refprop_TD('h', T0, rho0*(1-step), Y0))/(2*step*rho0);
%     dh_dY = (refprop_TD('h', T0, rho0, Y0*(1+step)) - refprop_TD('h', T0, rho0, Y0*(1-step)))/(2*step*Y0);
% 
%     lambda = refprop_TD('lambda', T0, rho0, Y0);
%     dlambda_dT = (refprop_TD('lambda', T0*(1+step), rho0, Y0) - refprop_TD('lambda', T0*(1-step), rho0, Y0))/(2*step*T0);
%     dlambda_drho = (refprop_TD('lambda', T0, rho0*(1+step), Y0) - refprop_TD('lambda', T0, rho0*(1-step), Y0))/(2*step*rho0);
%     dlambda_dY = (refprop_TD('lambda', T0, rho0, Y0*(1+step)) - refprop_TD('lambda', T0, rho0, Y0*(1-step)))/(2*step*Y0);
%     d2lambda_dTdT = (refprop_TD('lambda', T0*(1+step), rho0, Y0) - 2*D12 + refprop_TD('lambda', T0*(1-step), rho0, Y0))/(step*T0)^2;
%     d2lambda_dTdrho = (refprop_TD('lambda', T0*(1+step), rho0*(1+step), Y0)...
%         - refprop_TD('lambda', T0*(1+step), rho0*(1-step), Y0)...
%         - refprop_TD('lambda', T0*(1-step), rho0*(1+step), Y0)...
%         + refprop_TD('lambda', T0*(1-step), rho0*(1-step), Y0))/(4*step^2*T0*rho0);
%     d2lambda_dTdY = (refprop_TD('lambda', T0*(1+step), rho0, Y0*(1+step))...
%         - refprop_TD('lambda', T0*(1+step), rho0, Y0*(1-step))...
%         - refprop_TD('lambda', T0*(1-step), rho0, Y0*(1+step))...
%         + refprop_TD('lambda', T0*(1-step), rho0, Y0*(1-step)))/(4*step^2*T0*Y0);
%     d2lambda_drhodrho = (refprop_TD('lambda', T0, rho0*(1+step), Y0) - 2*D12 + refprop_TD('lambda', T0, rho0*(1-step), Y0))/(step*rho0)^2;
%     d2lambda_drhodY = (refprop_TD('lambda', T0, rho0*(1+step), Y0*(1+step))...
%         - refprop_TD('lambda', T0, rho0*(1+step), Y0*(1-step))...
%         - refprop_TD('lambda', T0, rho0*(1-step), Y0*(1+step))...
%         + refprop_TD('lambda', T0, rho0*(1-step), Y0*(1-step)))/(4*step^2*rho0*Y0);
%     d2lambda_dYdY = (refprop_TD('lambda', T0, rho0, Y0*(1+step)) - 2*D12 + refprop_TD('lambda', T0, rho0, Y0*(1-step)))/(step*Y0)^2;

    D12 = refprop_TD('d', T0, rho0, Y0, paramArray);
    [grad_d, err] = gradest(@(vec) refprop_TD('d', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end 
    [hess_d, err] = hessian(@(vec) refprop_TD('d', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    
    dD_dT = grad_d(1);
    dD_drho = grad_d(2);
    dD_dY = grad_d(3);
    d2D_dTdT = hess_d(1,1);
    d2D_dTdrho = hess_d(1,2);
    d2D_dTdY = hess_d(1,3);
    d2D_drhodrho = hess_d(2,2);
    d2D_drhodY = hess_d(2,3);
    d2D_dYdY = hess_d(3,3);
    
    
    mu = refprop_TD('mu', T0, rho0, Y0, paramArray);
    [grad_mu, err] = gradest(@(vec) refprop_TD('mu', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    [hess_mu, err] = hessian(@(vec) refprop_TD('mu', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    
    dmu_dT = grad_mu(1);
    dmu_drho = grad_mu(2);
    dmu_dY = grad_mu(3);
    d2mu_dTdT = hess_mu(1,1);
    d2mu_dTdrho = hess_mu(1,2);
    d2mu_dTdY = hess_mu(1,3);
    d2mu_drhodrho = hess_mu(2,2);
    d2mu_drhodY = hess_mu(2,3);
    d2mu_dYdY = hess_mu(3,3);
    
    P = refprop_TD('P', T0, rho0, Y0, paramArray);
    [grad_P, err] = gradest(@(vec) refprop_TD('P', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    [hess_P, err] = hessian(@(vec) refprop_TD('P', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    
    dp_dT = grad_P(1);
    dp_drho = grad_P(2);
    dp_dY = grad_P(3);
    d2p_dTdT = hess_P(1,1);
    d2p_dTdrho = hess_P(1,2);
    d2p_dTdY = hess_P(1,3);
    d2p_drhodrho = hess_P(2,2);
    d2p_drhodY = hess_P(2,3);
    d2p_dYdY = hess_P(3,3);

    h = refprop_TD('h', T0, rho0, Y0, paramArray);
    [grad_h, err] = gradest(@(vec) refprop_TD('h', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    
    dh_dT = grad_h(1);
    dh_drho = grad_h(2);
    dh_dY = grad_h(3);
    
    lambda = refprop_TD('lambda', T0, rho0, Y0, paramArray);
    [grad_lambda, err] = gradest(@(vec) refprop_TD('lambda', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    [hess_lambda, err] = hessian(@(vec) refprop_TD('lambda', vec(1), vec(2), vec(3), paramArray), [T0 rho0 Y0]);
    if max(err(:)) > 1e-6
        warning('Error = %.4g', max(err(:)));
    end    
    
    dlambda_dT = grad_lambda(1);
    dlambda_drho = grad_lambda(2);
    dlambda_dY = grad_lambda(3);
    d2lambda_dTdT = hess_lambda(1,1);
    d2lambda_dTdrho = hess_lambda(1,2);
    d2lambda_dTdY = hess_lambda(1,3);
    d2lambda_drhodrho = hess_lambda(2,2);
    d2lambda_drhodY = hess_lambda(2,3);
    d2lambda_dYdY = hess_lambda(3,3);
end