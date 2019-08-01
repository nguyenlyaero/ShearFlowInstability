function [cp_R] = PC_SAFT_DT_cp(rho_hat,T,x,S1,S2)
% USING FINITE-DIFFERENCING
% See: Lemmon 1999
% cp_R: constant pressure specific heat
%
% T: [K]
% rho_hat: molar density (kmol/m^3)
% x: molar fraction of first species
    function a_ig = getIdealHelmholtz_2S(tau)
            Ti = Tc/tau;
            %P = PC_SAFT_DT_Cubic(rho_hat,Ti,x,S1,S2);
            P = rho_hat*8.31*1000*Ti;
            h_ig_RT = IdealEnthalpy(Ti, x, S1, S2);
            s_ig_R = IdealEntropy(Ti, P, x, S1, S2);
            a_ig = h_ig_RT - 1 - s_ig_R;       
    end
    function a_ig = getIdealHelmholtz_1S(tau)
            Ti = Tc/tau;
            %P = PC_SAFT_DT_Cubic(rho_hat,Ti,x,S1);
            P = rho_hat*8.31*1000*Ti;
            h_ig_RT = IdealEnthalpy(Ti, x, S1);
            s_ig_R = IdealEntropy(Ti, P, x, S1);
            a_ig = h_ig_RT - 1 - s_ig_R;       
    end
    function a_r = getResHelmholtz_2S(tau, delta)
            Ti = Tc/tau;
            rho_hati = delta/Vc;
            [~,a_r,~,~] = PC_SAFT_DT_Cubic(rho_hati,Ti,x,S1,S2);   
    end
    function a_r = getResHelmholtz_1S(tau, delta)
            Ti = Tc/tau;
            rho_hati = delta/Vc;
            [~,a_r,~,~] = PC_SAFT_DT_Cubic(rho_hati,Ti,x,S1);    
    end
    if nargin == 5
        [Tc1, ~, ~, ~, ~ , ~, ~, ~, ~, Vc1] = getSpecies(S1);
        [Tc2, ~, ~, ~, ~ , ~, ~, ~, ~, Vc2] = getSpecies(S2);
        Tc = x*Tc1 + (1-x)*Tc2;
        Vc = x*Vc1 + (1-x)*Vc2; %[m^3/kmol]
        
        d2a0_dtau2 = derivest(@getIdealHelmholtz_2S, Tc/T, 'Vectorized', 'no', 'DerivativeOrder', 2);
        
        dar_ddelta = derivest(@(delta) getResHelmholtz_2S(Tc/T, delta), rho_hat*Vc, 'Vectorized', 'no');
        
        hess = hessian(@(var) getResHelmholtz_2S(var(1), var(2)), [Tc/T, rho_hat*Vc]);
        d2ar_dtau2 = hess(1,1);
        d2ar_ddelta2 = hess(2,2);
        d2ar_dtauddelta = hess(1,2);
        
        tau = Tc/T;
        delta = rho_hat*Vc;
        
        cv_R = -tau^2*(d2a0_dtau2 + d2ar_dtau2);
        cp_R = cv_R + (1 + delta*dar_ddelta - delta*tau*d2ar_dtauddelta)^2/ ...
            (1 + 2*delta*dar_ddelta + delta^2*d2ar_ddelta2);
    elseif nargin == 4
        [Tc, ~, ~, ~, ~ , ~, ~, ~, ~, Vc] = getSpecies(S1);
        
        d2a0_dtau2 = derivest(@getIdealHelmholtz_1S, Tc/T, 'Vectorized', 'no', 'DerivativeOrder', 2);
        
        dar_ddelta = derivest(@(delta) getResHelmholtz_1S(Tc/T, delta), rho_hat*Vc, 'Vectorized', 'no');
        
        hess = hessian(@(var) getResHelmholtz_1S(var(1), var(2)), [Tc/T, rho_hat*Vc]);
        d2ar_dtau2 = hess(1,1);
        d2ar_ddelta2 = hess(2,2);
        d2ar_dtauddelta = hess(1,2);
        
        tau = Tc/T;
        delta = rho_hat*Vc;
        
        cv_R = -tau^2*(d2a0_dtau2 + d2ar_dtau2);
        cp_R = cv_R + (1 + delta*dar_ddelta - delta*tau*d2ar_dtauddelta)^2/ ...
            (1 + 2*delta*dar_ddelta + delta^2*d2ar_ddelta2);
    end
end

