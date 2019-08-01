function [Viscosity, Lambda] = Chung_PT(P, T, x, S1, S2)
    % Get Species Data
    if nargin==4
        [Tc1, ~, acentric1, M1, ~ , ~, ~, ~, dipole1, Vc1] = getSpecies(S1);        
        
        acentric_vec = [acentric1; acentric1];
        W_vec = [M1; M1];
        epsilon_k_vec = [Tc1/1.2593; Tc1/1.2593];
        sigma_vec = [0.809*(1000*Vc1)^(1/3); 0.809*(1000*Vc1)^(1/3)];
        dipole_vec = [dipole1; dipole1];
        x_vec = [0.5; 0.5];
        
        rho_hat = PC_SAFT_PT_Cubic(P,T,x,S1,0.5);
        [~, cp_ig_R] = IdealEnthalpy(T, x, S1);
    elseif nargin==5
        [Tc1, ~, acentric1, M1, ~ , ~, ~, ~, dipole1, Vc1] = getSpecies(S1);
        [Tc2, ~, acentric2, M2, ~ , ~, ~, ~, dipole2, Vc2] = getSpecies(S2);
        
        acentric_vec = [acentric1; acentric2];
        W_vec = [M1; M2];
        epsilon_k_vec = [Tc1/1.2593; Tc2/1.2593];
        sigma_vec = [0.809*(1000*Vc1)^(1/3); 0.809*(1000*Vc2)^(1/3)];
        dipole_vec = [dipole1; dipole2];
        x_vec = [x, 1-x];    
        
        rho_hat = PC_SAFT_PT_Cubic(P,T,x,S1,S2,0.5);
        [~, cp_ig_R] = IdealEnthalpy(T, x, S1, S2);
    else
        assert(0, 'Wrong number of inputs');
    end
    
    % Mixture average
    W_ij = zeros(2,2);
    omega_ij = zeros(2,2);
    sigma_ij = zeros(2,2);
    epsilon_k_ij = zeros(2,2);
    for i=1:2
        for j=1:2
            W_ij(i,j) = 2*W_vec(i)*W_vec(j)/(W_vec(i) + W_vec(j));
            omega_ij(i,j) = 1/2*(acentric_vec(i) + acentric_vec(j));
            sigma_ij(i,j) = sqrt(sigma_vec(i)*sigma_vec(j));
            epsilon_k_ij(i,j) = sqrt(epsilon_k_vec(i)*epsilon_k_vec(j));
        end
    end
    W_m = 0;
    omega_m = 0;
    dipole_m = 0;
    sigma_m = 0;
    epsilon_k_m = 0;
    for i=1:2
        for j=1:2
            W_m = W_m + x_vec(i)*x_vec(j)*epsilon_k_ij(i,j)*sigma_ij(i,j)^2*sqrt(W_ij(i,j));
            omega_m = omega_m + x_vec(i)*x_vec(j)*omega_ij(i,j)*sigma_ij(i,j)^3;
            dipole_m = dipole_m + x_vec(i)*x_vec(j)*dipole_vec(i)^2*dipole_vec(j)^2/sigma_ij(i,j)^3;
            sigma_m = sigma_m + x_vec(i)*x_vec(j)*sigma_ij(i,j)^3;
            epsilon_k_m = epsilon_k_m + x_vec(i)*x_vec(j)*epsilon_k_ij(i,j)*sigma_ij(i,j)^3;
        end
    end
    sigma_m = sigma_m^(1/3);
    epsilon_k_m = epsilon_k_m/sigma_m^3;
    W_m = (W_m/(epsilon_k_m*sigma_m^2))^2;
    omega_m = omega_m/sigma_m^3;
    dipole_m = (sigma_m^3*dipole_m)^(1/4);
    

    Vc_m = (sigma_m/0.809)^3;
    Tc_m = 1.2593*epsilon_k_m;
    dipole_r = 131.3*dipole_m/sqrt(Vc_m*Tc_m);
    T_star_m = T/epsilon_k_m;
    
    % Neufeld
    Neufeld = 1.16145*T_star_m^(-0.14874) + 0.52487*exp(-0.7732*T_star_m) ...
        + 2.16178*exp(-2.43787*T_star_m); % Eq (9-4-3) in Poling
    kappa_m = 0;
    Fc_m = 1.0 - 0.2756*omega_m + 0.059035*dipole_r^4 + kappa_m;
    
    % Low-Pressure Viscosity
    mu = 40.785*Fc_m*sqrt(W_m*T)/(Neufeld*Vc_m^(2/3))/10.0e6;
    
    % High-Pressure Viscosity
    a1 = [6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900, 24.275, 0.7972, -0.2382, 0.06863];
    b1 = [50.412, -1.154e-3, 254.209, 38.096, 7.630, -12.537, 3.450, 1.1170, 0.0677, 0.3479];
    c1 = [-51.680, -6.257e-3, -168.480, -8.464, -14.354, 4.985, -11.291, 0.01235, -0.8163, 0.5926];
    d1 = [1189.000, 0.03728, 3898.0, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025, -0.727];
    E = zeros(10,1);
    for i=1:10
        E(i) = a1(i) + b1(i)*omega_m + c1(i)*dipole_r^4 + d1(i)*kappa_m;
    end
    
    
    y = rho_hat/1000*Vc_m/6;
    G1 = (1.0 - 0.5*y)/(1.0 - y)^3;
    G2 = (E(1)/y*(1 - exp(-E(4)*y)) + E(2)*G1*exp(E(5)*y) ...
        + E(3)*G1)/(E(1)*E(4) + E(2) + E(3));
    eta_star_star = E(7)*y^2*G2*exp(E(8) + E(9)/T_star_m + E(10)/T_star_m^2);
    eta_star = T_star_m^0.5/Neufeld*(Fc_m*(1.0/G2 + E(6)*y)) + eta_star_star;
    
    Viscosity = (36.644*eta_star*sqrt(W_m*Tc_m)/Vc_m^(2/3))/10.0e6;
    
    % Thermal Conductivity
    a2 = [2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089];
    b2 = [0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11];
    c2 = [-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -54.217];
    d2 = [121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81];
    B = zeros(10,1);
    for i=1:7
        B(i) = a2(i) + b2(i)*omega_m + c2(i)*dipole_r^4 + d2(i)*kappa_m;
    end
    G2 = ((B(1)/y)*(1 - exp(-B(4)*y)) + B(2)*G1*exp(B(5)*y) + B(3)*G1) ...
        /(B(1)*B(4) + B(2) + B(3));
    
    cv_ig_R = cp_ig_R - 1;
    alpha = cv_ig_R - 3/2;
    beta = 0.7862 - 0.7109*omega_m + 1.3168*omega_m^2;
    Zr = 2.0 + 10.5*(T/Tc_m)^2;
    q = 0.003586*sqrt(Tc_m/(W_m/1000))/Vc_m^(2/3);
    psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665 * Zr) ...
            /(0.6366 + beta*Zr + 1.061*alpha*beta));
    Lambda = 31.2*mu*psi/(W_m/1000)*(1/G2 + B(6)*y) + q*B(7)*y^2*sqrt(T/Tc_m)*G2;
              
end

