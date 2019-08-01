function [x_vec, y_vec, L_F] = PC_SAFT_PT_FLASH(P,T,z,S1,S2)
% a: reduced Helholtz free energy (a/RT)
% rho_hat: molar density (kmol/m^3)
%
% T: [K]
% P: Pressure (Pa)
% z: molar fraction of first species
% str: species (e.g. 'H2&O2')
% eta0: initial guess, 1e-10 for vapor, 0.5 for liquid

    MaxIter = 1000;
    RelTol = 1e-8;

    N_Av=6.022e23;

    %% Get Species Data
    [Tc1, Pc1, acentric1, M1, m1 , epsilon_k1, sigma1] = getSpecies(S1);
    [Tc2, Pc2, acentric2, M2, m2 , epsilon_k2, sigma2] = getSpecies(S2);
    m_vec = [m1; m2];
    epsilon_k_vec = [epsilon_k1; epsilon_k2];
    sigma_vec = [sigma1; sigma2];
    x_vec = [x; 1-x];
    k_bin = getBinaryInteraction(S1, S2);

    Nspecies = length(m_vec);
    

    %%
    K_vec = zeros(Nspecies,1);
    x_vec = z_vec;
    y_vec = z_vec;
    L_F = 0.5;
    function [out, dout] = OBJ(L_F)
        step = 1e-30;
        out = sum(z_vec.*(1-K_vec)./(K_vec+L_F*(1-K_vec)));
        dout = imag(sum(z_vec.*(1-K_vec)./(K_vec+L_F*complex(1,step)*(1-K_vec))))/step/L_F;
    end
    % First Iteration
    % Estimate initial K value
    for i=1:Nspecies
        K_vec(i) = Pc_vec(i)*10^(7/3*(1+omega_vec(i))*(1-Tc_vec(i)/T))/P;
    end
    % Solve for next L_F
    next_L_F = L_F;
    [out, dout] = OBJ(next_L_F);
    iter1 = 0;
    while iter1 < MaxIter && abs(out) > RelTol
        iter1 = iter1 + 1;
        next_L_F = next_L_F - out/dout;
        [out, dout] = OBJ(next_L_F);
    end
    assert(abs(out)<RelTol);
    % Next x and y
    next_x_vec = z_vec./(K_vec+next_L_F*(1-K_vec));
    next_y_vec = next_x_vec.*K_vec;
    % Difference
    diff_L_F = abs((next_L_F - L_F)/L_F);
    diff_x = norm(next_x_vec - x_vec)/norm(x_vec);
    diff_y = norm(next_y_vec - y_vec)/norm(y_vec);
    
    iter2 = 0;
    while iter2 < MaxIter && (diff_L_F > RelTol || diff_x > RelTol || diff_y > RelTol)
        iter2 = iter2 + 1;
        L_F = next_L_F;
        x_vec = next_x_vec;
        y_vec = next_y_vec;
        % Find K value
        phi_vec_VAPOR = zeros(2,2);
        fug_VAPOR = zeros(2,1);
        phi_vec_LIQUID = zeros(2,2);
        fug_LIQUID = zeros(2,1);
            [~,~,log_phi_vec_VAPOR] = PC_SAFT_PT_Cubic(P,T,y_vec(1),S1,S2,1e-10);
            phi_vec_VAPOR(:,1) = exp(log_phi_vec_VAPOR);
            [~,~,log_phi_vec_VAPOR] = PC_SAFT_PT_Cubic(P,T,y_vec(1),S1,S2,0.5);
            phi_vec_VAPOR(:,2) = exp(log_phi_vec_VAPOR);
            fug_VAPOR = sum(y_vec.*phi_vec_VAPOR);
            [~,imin_VAPOR] = min(fug_VAPOR);
            
            [~,~,log_phi_vec_LIQUID] = PC_SAFT_PT_Cubic(P,T,x_vec(1),S1,S2,1e-10);
            phi_vec_LIQUID(:,1) = exp(log_phi_vec_LIQUID);
            [~,~,log_phi_vec_LIQUID] = PC_SAFT_PT_Cubic(P,T,x_vec(1),S1,S2,0.5);
            phi_vec_LIQUID(:,2) = exp(log_phi_vec_LIQUID);
            fug_LIQUID = sum(x_vec.*phi_vec_LIQUID);
            [~,imin_LIQUID] = min(fug_LIQUID);

        K_vec = phi_vec_LIQUID(:,imin_LIQUID)./phi_vec_VAPOR(:,imin_VAPOR);
        % Solve for next L_F
        next_L_F = L_F;
        [out, dout] = OBJ(next_L_F);
        iter1 = 0;
        while iter1 < MaxIter && abs(out) > RelTol
            iter1 = iter1 + 1;
            next_L_F = next_L_F - out/dout;
            [out, dout] = OBJ(next_L_F); 
        end
        assert(abs(out)<RelTol);
        % Next x and y
        next_x_vec = z_vec./(K_vec+next_L_F*(1-K_vec));
        next_y_vec = next_x_vec.*K_vec;
        % Difference
        diff_L_F = abs((next_L_F - L_F)/L_F);
        diff_x = norm(next_x_vec - x_vec)/norm(x_vec);
        diff_y = norm(next_y_vec - y_vec)/norm(y_vec);
    end
end

