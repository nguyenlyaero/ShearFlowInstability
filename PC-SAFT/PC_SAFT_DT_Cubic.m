function [P,a_h,log_phi_vec,h_hat] = PC_SAFT_DT_Cubic(rho_hat,T,x,S1,S2)
% a: reduced Helholtz free energy (a/RT)
% P: Pressure (Pa)
% log_phi_vec: log of fugacity coefficient for each species
%
% T: [K]
% rho_hat: molar density (kmol/m^3)
% x: molar fraction of first species
% str: species (e.g. 'H2&O2')
% h_hat: molar enthalpy J/mol 
    R = 8.314462618; % J/mol K
    N_Av=6.022e23;
    rho=rho_hat.*N_Av*10^3/(10^10)^3; % number molecular density (1/Angstrom)^3

    %% Get Species Data
    if nargin == 4
        [Tc, Pc, acentric, M, m , epsilon_k, sigma] = getSpecies(S1);
        m_vec = m;
        epsilon_k_vec = epsilon_k;
        sigma_vec = sigma;
        x_vec = 1;
        k_bin = 0;
        
        h_ig_RT = IdealEnthalpy(T, x, S1);
    else % nargin == 5
        [Tc1, Pc1, acentric1, M1, m1 , epsilon_k1, sigma1] = getSpecies(S1);
        [Tc2, Pc2, acentric2, M2, m2 , epsilon_k2, sigma2] = getSpecies(S2);
        m_vec = [m1; m2];
        epsilon_k_vec = [epsilon_k1; epsilon_k2];
        sigma_vec = [sigma1; sigma2];
        x_vec = [x; 1-x];
        k_bin = getBinaryInteraction(S1, S2);
        
        h_ig_RT = IdealEnthalpy(T, x, S1, S2);
    end
    
    Nspecies = length(m_vec);
    
    %% Pressure
    m_bar = sum(x_vec.*m_vec);

    d_vec = sigma_vec.*(1-0.12.*exp(-3*epsilon_k_vec/T)); % Temp. Dependent segment diameter

    Xis=zeros(4,1); % Reduced densities
    for i=1:4
        Xis(i) = pi/6*rho*sum(x_vec.*m_vec.*d_vec.^(i-1));
    end

    % Radial Distribution Function
    g_hs_vec = 1/(1-Xis(4)) + d_vec/2*3*Xis(3)/(1-Xis(4))^2 + d_vec.^2/4*2*Xis(3)^2/(1-Xis(4))^3;

    eta=Xis(4);
    % Model Constants
    alpha=zeros(3,7);
    alpha(1,1)=0.9105631445; alpha(2,1)=-0.3084016918; alpha(3,1)=-0.0906148351;
    alpha(1,2)=0.6361281449; alpha(2,2)=0.1860531159; alpha(3,2)=0.4527842806;
    alpha(1,3)=2.6861347891; alpha(2,3)=-2.5030047259; alpha(3,3)=0.5962700728;
    alpha(1,4)=-26.547362491; alpha(2,4)=21.419793629; alpha(3,4)=-1.7241829131;
    alpha(1,5)=97.759208784; alpha(2,5)=-65.255885330; alpha(3,5)=-4.1302112531;
    alpha(1,6)=-159.59154087; alpha(2,6)=83.318680481; alpha(3,6)=13.776631870;
    alpha(1,7)=91.297774084; alpha(2,7)=-33.746922930; alpha(3,7)=-8.6728470368;
    beta=zeros(3,7);
    beta(1,1)=0.7240946941; beta(2,1)=-0.5755498075; beta(3,1)=0.0976883116;
    beta(1,2)=2.2382791861; beta(2,2)=0.6995095521; beta(3,2)=-0.2557574982;
    beta(1,3)=-4.0025849485; beta(2,3)=3.8925673390; beta(3,3)=-9.1558561530;
    beta(1,4)=-21.003576815; beta(2,4)=-17.215471648; beta(3,4)=20.642075974;
    beta(1,5)=26.855641363; beta(2,5)=192.67226447; beta(3,5)=-38.804430052;
    beta(1,6)=206.55133841; beta(2,6)=-161.82646165; beta(3,6)=93.626774077;
    beta(1,7)=-355.60235612; beta(2,7)=-165.20769346; beta(3,7)=-29.666905585;
    a=zeros(7,1);
    b=zeros(7,1);
    for i=1:7
        a(i) = alpha(1,i) + (m_bar-1)/m_bar*alpha(2,i) + (m_bar-1)/m_bar*(m_bar-2)/m_bar*alpha(3,i);
        b(i) = beta(1,i) + (m_bar-1)/m_bar*beta(2,i) + (m_bar-1)/m_bar*(m_bar-2)/m_bar*beta(3,i);
    end

    % Integral Models
    I1 = 0;
    I2 = 0;
    for i=1:7
        I1 = I1+a(i)*eta^(i-1);
        I2 = I2+b(i)*eta^(i-1);
    end

    % Compressibility Term
    C_1 = (1 + m_bar*(8*eta-2*eta^2)/(1-eta)^4 ...
        + (1-m_bar)*(20*eta-27*eta^2+12*eta^3-2*eta^4)/((1-eta)*(2-eta))^2)^-1;
    % Mixing Rules
    sigma_mat = zeros(Nspecies,Nspecies);
    epsilon_k_mat = zeros(Nspecies,Nspecies);
    for i=1:Nspecies
        for j=1:Nspecies
            sigma_mat(i,j) = 1/2*(sigma_vec(i)+sigma_vec(j));
            if i~=j
                epsilon_k_mat(i,j) = sqrt(epsilon_k_vec(i)*epsilon_k_vec(j))*(1-k_bin);
            else
                epsilon_k_mat(i,j) = epsilon_k_vec(i);
            end
        end
    end

    m2_epsilon_k_sigma_3_bar = 0;
    m2_epsilon_k_2_sigma_3_bar = 0;
    for i=1:Nspecies
        for j=1:Nspecies
            m2_epsilon_k_sigma_3_bar = m2_epsilon_k_sigma_3_bar + ...
                x_vec(i)*x_vec(j)*m_vec(i)*m_vec(j)*epsilon_k_mat(i,j)*sigma_mat(i,j)^3;
            m2_epsilon_k_2_sigma_3_bar = m2_epsilon_k_2_sigma_3_bar + ...
                x_vec(i)*x_vec(j)*m_vec(i)*m_vec(j)*epsilon_k_mat(i,j)^2*sigma_mat(i,j)^3;
        end
    end

    % Dispersion energy contribution
    a_disp = -2*pi*rho*m2_epsilon_k_sigma_3_bar/T*I1 ...
        - pi*rho*m_bar*C_1*m2_epsilon_k_2_sigma_3_bar/T^2*I2;

    

    % Hard-Chain Compressibility Contribution
    Z_hs = Xis(4)/(1-Xis(4)) + 3*Xis(2)*Xis(3)/Xis(1)/(1-Xis(4))^2 ...
        + (3*Xis(3)^3 - Xis(4)*Xis(3)^3)/Xis(1)/(1-Xis(4))^3;
    B_vec = Xis(4)/(1-Xis(4))^2 + d_vec/2*(3*Xis(3)/(1-Xis(4))^2 + 6*Xis(3)*Xis(4)/(1-Xis(4))^3) ...
        + d_vec.^2/4*(4*Xis(3)^2/(1-Xis(4))^3 + 6*Xis(3)^2*Xis(4)/(1-Xis(4))^4);
    Z_hc = m_bar*Z_hs-sum(x_vec.*(m_vec-1)./g_hs_vec.*B_vec);

    % Dispersion Compressibility Contribution
    I1Diff=0; I2Diff=0;
    for i=1:7
        I1Diff = I1Diff+a(i).*i.*eta.^(i-1);
        I2Diff = I2Diff+b(i).*i.*eta.^(i-1);
    end
    C_2 = -C_1^2*(m_bar*(-4*eta^2+20*eta+8)/(1-eta)^5 ...
        + (1-m_bar)*(2*eta^3+12*eta^2-48*eta+40)/((1-eta)*(2-eta))^3);

    Z_disp = -2*pi*rho*I1Diff*m2_epsilon_k_sigma_3_bar/T ...
        - pi*rho*m_bar*(C_1*I2Diff+C_2*eta*I2)*m2_epsilon_k_2_sigma_3_bar/T^2;

    % Compressibility Factor
    Z = 1 + Z_hc + Z_disp;

    % Pressure
    k = 1.380649e-23; % Boltzmann Constant (J/K)
    P = Z.*k.*T.*rho.*(10^10)^3;

    %% Residual Helholtz Energy
    % Hard chain energy contribution
    a_hs = 1/Xis(1) * (3*Xis(2)*Xis(3)/(1-Xis(4)) + Xis(3)^3./(Xis(4)*(1-Xis(4))^2) ...
        + (Xis(3)^3/Xis(4)^2-Xis(1))*log(1-Xis(4)));
    a_hc = m_bar*a_hs-sum(x_vec.*(m_vec-1).*log(g_hs_vec));
    % Residual Helholtz Energy
    a_h = a_hc+a_disp;
    
    %% Fugacity Coefficient
    dXi_n_d_x_k = zeros(4,Nspecies);
    for i=1:4
        for j=1:Nspecies
            dXi_n_d_x_k(i,j) = pi/6*rho*m_vec(j)*d_vec(j)^(i-1);
        end
    end
    
    % Hard-Chain Contribution
    dg_i_dx_k = zeros(Nspecies,Nspecies);
    for i=1:Nspecies        
        for j=1:Nspecies
            dg_i_dx_k(i,j) = dXi_n_d_x_k(4,j)/(1-Xis(4))^2 ...
                + d_vec(i)/2*(3*dXi_n_d_x_k(3,j)/(1-Xis(4))^2 + 6*Xis(3)*dXi_n_d_x_k(4,j)/(1-Xis(4)^3)) ...
                + d_vec(i)^2/4*(4*Xis(3)*dXi_n_d_x_k(3,j)/(1-Xis(4))^3 + 6*Xis(3)^2*dXi_n_d_x_k(4,j)/(1-Xis(4))^4);
        end
    end

    dahs_dx_k = zeros(Nspecies,1);
    for n=1:Nspecies
        dahs_dx_k(n) = -dXi_n_d_x_k(1,n)/Xis(1)*a_hs + 1/Xis(1)*( ...
            3*(dXi_n_d_x_k(2,n)*Xis(3) + Xis(2)*dXi_n_d_x_k(3,n))/(1-Xis(4)) ...
            + 3*Xis(2)*Xis(3)*dXi_n_d_x_k(4,n)/(1-Xis(4))^2 + 3*Xis(3)^2*dXi_n_d_x_k(3,n)/Xis(4)/(1-Xis(4))^2 ...
            + Xis(3)^3*dXi_n_d_x_k(4,n)*(3*Xis(4)-1)/Xis(4)^2/(1-Xis(4))^3 ...
            + ((3*Xis(3)^2*dXi_n_d_x_k(3,n)*Xis(4)-2*Xis(3)^3*dXi_n_d_x_k(4,n))/Xis(4)^3-dXi_n_d_x_k(1,n))*log(1-Xis(4)) ...
            + (Xis(1)-Xis(3)^3/Xis(4)^2)*dXi_n_d_x_k(4,n)/(1-Xis(4)));
    end
    
    da_hc_dx_k = zeros(Nspecies,1);
    for n=1:Nspecies
        da_hc_dx_k(n) = m_vec(n)*a_hs + m_bar*dahs_dx_k(n) ...
            - sum(x_vec.*(m_vec-1)./g_hs_vec.*dg_i_dx_k(:,n));
    end
    
    % Dispersion Contribution
    da_i_dx_k = zeros(7,Nspecies);
    db_i_dx_k = zeros(7,Nspecies);
    for i=1:7
        for n=1:Nspecies
            da_i_dx_k(i,n) = m_vec(n)/m_bar^2*alpha(2,i) + m_vec(n)/m_bar^2*(3-4/m_bar)*alpha(3,i);
            db_i_dx_k(i,n) = m_vec(n)/m_bar^2*beta(2,i) + m_vec(n)/m_bar^2*(3-4/m_bar)*beta(3,i);
        end
    end
    
    dI1_dx_k = zeros(Nspecies,1);
    dI2_dx_k = zeros(Nspecies,1);
    for n=1:Nspecies
        for i=1:7
            dI1_dx_k(n) = dI1_dx_k(n) + (i-1)*a(i)*dXi_n_d_x_k(4,n)*eta^(i-2) + da_i_dx_k(i,n)*eta^(i-1);
            dI2_dx_k(n) = dI2_dx_k(n) + (i-1)*b(i)*dXi_n_d_x_k(4,n)*eta^(i-2) + db_i_dx_k(i,n)*eta^(i-1);
        end
    end
    
    dC1_dx_k = zeros(Nspecies,1);
    for n=1:Nspecies
        dC1_dx_k(n) = C_2*dXi_n_d_x_k(4,n) ...
            - C_1^2*(m_vec(n)*(8*eta-2*eta^2)/(1-eta)^4 ...
            - m_vec(n)*(20*eta-27*eta^2+12*eta^3-2*eta^4)/((1-eta)*(2-eta))^2);
    end
    
    dm2_epsilon_k_sigma_3_bar_dx_k = zeros(Nspecies,1);
    dm2_epsilon_k2_sigma_3_bar_dx_k = zeros(Nspecies,1);
    for i=1:Nspecies
        for j = 1:Nspecies
            dm2_epsilon_k_sigma_3_bar_dx_k(i) = ...
                dm2_epsilon_k_sigma_3_bar_dx_k(i) + 2*m_vec(i)*x_vec(j)*m_vec(j)*epsilon_k_mat(i,j)*sigma_mat(i,j)^3;
            dm2_epsilon_k2_sigma_3_bar_dx_k(i) = ...
                dm2_epsilon_k2_sigma_3_bar_dx_k(i) + 2*m_vec(i)*x_vec(j)*m_vec(j)*epsilon_k_mat(i,j)^2*sigma_mat(i,j)^3;
        end        
    end
    
    da_disp_dx_k = zeros(Nspecies,1);
    for n=1:Nspecies
        da_disp_dx_k(n) = -2*pi*rho*(dI1_dx_k(n)*m2_epsilon_k_sigma_3_bar/T + I1*dm2_epsilon_k_sigma_3_bar_dx_k(n)/T) ...
            - pi*rho*((m_vec(n)*C_1*I2 + m_bar*dC1_dx_k(n)*I2 + m_bar*C_1*dI2_dx_k(n))*m2_epsilon_k_2_sigma_3_bar/T^2 ...
            + m_bar*C_1*I2*dm2_epsilon_k2_sigma_3_bar_dx_k(n)/T^2);
    end
    
    da_h_dx_k = da_hc_dx_k + da_disp_dx_k;
    
    mu_kT_vec = zeros(Nspecies,1);
    for n=1:Nspecies
        mu_kT_vec(n) = a_h + (Z-1) + da_h_dx_k(n) - sum(x_vec.*da_h_dx_k);
    end
    
    log_phi_vec = zeros(Nspecies,1);
    for n=1:Nspecies
        log_phi_vec(n) = mu_kT_vec(n)-log(Z);
    end
    
    %% Enthalpy
    dd_dT_vec = sigma_vec.*(3*epsilon_k_vec/T^2).*(-0.12*exp(-3*epsilon_k_vec/T));
    dXis_dT = zeros(4,1);
    for i=2:4
        dXis_dT(i) = pi/6*rho*sum(x_vec.*m_vec.*(i-1).*dd_dT_vec.*d_vec.^(i-2));
    end
    
    
    da_hs_dT = 1/Xis(1)*(3*dXis_dT(2)*Xis(3)/(1-Xis(4)) + 3*Xis(2)*Xis(3)*dXis_dT(4)/(1-Xis(4))^2 ...
        + 3*Xis(3)^2*dXis_dT(3)/Xis(4)/(1-Xis(4))^2 + Xis(3)^3*dXis_dT(4)*(3*Xis(4)-1)/Xis(4)^2/(1-Xis(4))^3 ...
        + (3*Xis(3)^2*dXis_dT(3)*Xis(4) - 2*Xis(3)^3*dXis_dT(4))/Xis(4)^3*log(1-Xis(4)) + (Xis(1)-Xis(3)^3/Xis(4)^2)*dXis_dT(4)/(1-Xis(4)));
    dg_hs_vec_dT = dXis_dT(4)/(1-Xis(4))^2 + dd_dT_vec/2*3*Xis(3)/(1-Xis(4))^2 ...
        + d_vec/2*(3*dXis_dT(3)/(1-Xis(4))^2 + 6*Xis(3)*dXis_dT(4)/(1-Xis(4))^3) ...
        + 1/2*d_vec.*dd_dT_vec*2*Xis(3)^2/(1-Xis(4))^3 + (1/2*d_vec).^2*(4*Xis(3)*dXis_dT(3)/(1-Xis(4))^3 + 6*Xis(3)^2*dXis_dT(4)/(1-Xis(4))^4);
    da_hc_dT = m_bar*da_hs_dT - sum(x_vec.*(m_vec-1)./g_hs_vec.*dg_hs_vec_dT);
    
    dI1_dT = 0;
    dI2_dT = 0;
    for i=1:7
        dI1_dT = dI1_dT + a(i)*(i-1)*dXis_dT(4)*eta^(i-2);
        dI2_dT = dI2_dT + b(i)*(i-1)*dXis_dT(4)*eta^(i-2);
    end
    dC1_dT = dXis_dT(4)*C_2;
    da_disp_dT = -2*pi*rho*(dI1_dT-I1/T)*m2_epsilon_k_sigma_3_bar/T ...
        - pi*rho*m_bar*(dC1_dT*I2 + C_1*dI2_dT - 2*C_1*I2/T)*m2_epsilon_k_2_sigma_3_bar/T^2; 
    dah_dT = da_hc_dT + da_disp_dT;
    h_res_RT = -T*dah_dT + (Z-1);
    
    h_hat_RT = h_ig_RT + h_res_RT;
    h_hat = h_hat_RT*R*T;
end

