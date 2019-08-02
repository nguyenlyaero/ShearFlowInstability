%% Load data
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
N = 800;
alpha = 0.58;
[xi_vec, D0, D1, D2, ~] = Dmat(N);


%% Plot Eigenspectrum
figure;
plot(real(omega_vec), imag(omega_vec),'o');
xlim([0.4 1.2]);
ylim([-0.7 0.1]);
xlabel('\omega_r');
ylabel('\omega_i');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g', Rey, N, alpha));

%% Calculate Profile
ii = 750;
eigen_vec_i = V(:,ii);

u_profile = zeros(N,1);
v_profile = zeros(N,1);
T_profile = zeros(N,1);
rho_profile = zeros(N,1);
Y_profile = zeros(N,1);
for j=1:N
    u_profile = u_profile + eigen_vec_i(j)*D0(:,j);
    v_profile = v_profile + eigen_vec_i(N+j)*D0(:,j);
    T_profile = T_profile + eigen_vec_i(2*N+j)*D0(:,j);
    rho_profile = rho_profile + eigen_vec_i(3*N+j)*D0(:,j);  
    Y_profile = Y_profile + eigen_vec_i(4*N+j)*D0(:,j);
end

p_profile = dp_dT_vec.*T_profile + dp_drho_vec.*rho_profile;
h_profile = dh_dT_vec.*T_profile + dh_drho_vec.*rho_profile;
mu_profile = dmu_dT_vec.*T_profile + dmu_drho_vec.*rho_profile;
lambda_profile = dlambda_dT_vec.*T_profile + dlambda_drho_vec.*rho_profile;
D_profile = dD_dT_vec.*T_profile + dD_drho_vec.*rho_profile;

figure;
hold on;
plot(xi_vec, abs(u_profile));
plot(xi_vec, real(u_profile), 'k--');
plot(xi_vec, imag(u_profile), 'k:');
hold off;
xlabel('\xi');
ylabel('|u_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

%% Eigenmode in Physical Space
N_zeta = 1000;
zeta_vec = linspace(0,4*pi/alpha,N_zeta);

u_profile_mat = zeros(N,N_zeta);
v_profile_mat = zeros(N,N_zeta);
h_profile_mat = zeros(N,N_zeta);
rho_profile_mat = zeros(N,N_zeta);
Y_profile_mat = zeros(N,N_zeta);

p_profile_mat = zeros(N,N_zeta);
T_profile_mat = zeros(N,N_zeta);
mu_profile_mat = zeros(N,N_zeta);
lambda_profile_mat = zeros(N,N_zeta);
D_profile_mat = zeros(N,N_zeta);

for j=1:N_zeta
    u_profile_mat(:,j) = real(u_profile*exp(1i*alpha*zeta_vec(j)));
    v_profile_mat(:,j) = real(v_profile*exp(1i*alpha*zeta_vec(j)));
    h_profile_mat(:,j) = real(h_profile*exp(1i*alpha*zeta_vec(j)));
    rho_profile_mat(:,j) = real(rho_profile*exp(1i*alpha*zeta_vec(j)));
    Y_profile_mat(:,j) = real(Y_profile*exp(1i*alpha*zeta_vec(j)));

    p_profile_mat(:,j) = real(p_profile*exp(1i*alpha*zeta_vec(j)));
    T_profile_mat(:,j) = real(T_profile*exp(1i*alpha*zeta_vec(j)));
    mu_profile_mat(:,j) = real(mu_profile*exp(1i*alpha*zeta_vec(j)));
    lambda_profile_mat(:,j) = real(lambda_profile*exp(1i*alpha*zeta_vec(j)));
    D_profile_mat(:,j) = real(D_profile*exp(1i*alpha*zeta_vec(j)));
end

[zeta_mat, xi_mat] = meshgrid(zeta_vec, xi_vec);
y_mat = xi_mat./(1-xi_mat.^2).^(1/2);

%% Basic Variables Perturbation
figure;

subplot(5,3,1);
plot1 = surface(alpha*zeta_mat, y_mat, u_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('u1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,2);
plot(u0__x_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('u_0');
subplot(5,3,3)
hold on;
plot(abs(u_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(u_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(u_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
xlabel('|u_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,4);
plot1 = surface(alpha*zeta_mat, y_mat, v_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('v1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,5);
plot(zeros(N,1), xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('v_0');
subplot(5,3,6)
hold on;
plot(abs(v_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(v_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(v_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
xlabel('|v_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,7);
plot1 = surface(alpha*zeta_mat, y_mat, h_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('h1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,8);
plot(h_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('h_0');
subplot(5,3,9)
hold on;
plot(abs(h_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(h_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(h_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
xlabel('|h_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,10);
plot1 = surface(alpha*zeta_mat, y_mat, rho_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('rho1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,11);
plot(rho0_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('rho_0');
subplot(5,3,12)
hold on;
plot(abs(rho_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(rho_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(rho_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
xlabel('|rho_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,13);
plot1 = surface(alpha*zeta_mat, y_mat, Y_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('Y1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,14);
plot(Y0_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('Y_0');
subplot(5,3,15)
hold on;
plot(abs(Y_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(Y_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(Y_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
xlabel('|Y_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

%% Derived Propertires Perturbation

figure;

subplot(5,3,1);
plot1 = surface(alpha*zeta_mat, y_mat, p_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('p1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,2);
plot(ones(N,1), xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('p_0');
subplot(5,3,3)
hold on;
plot(abs(p_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(p_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(p_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
ylabel('|p_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,4);
plot1 = surface(alpha*zeta_mat, y_mat, T_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('T1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,5);
plot(T0_vec*T_ref/T_base, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('T_0');
subplot(5,3,6)
hold on;
plot(abs(T_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(T_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(T_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
ylabel('|T_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,7);
plot1 = surface(alpha*zeta_mat, y_mat, mu_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('mu1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,8);
plot(mu_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('mu_0');
subplot(5,3,9)
hold on;
plot(abs(mu_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(mu_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(mu_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
ylabel('|mu_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,10);
plot1 = surface(alpha*zeta_mat, y_mat, lambda_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('lambda1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,11);
plot(lambda_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('lambda_0');
subplot(5,3,12)
hold on;
plot(abs(lambda_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(lambda_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(lambda_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
ylabel('|lambda_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));

subplot(5,3,13);
plot1 = surface(alpha*zeta_mat, y_mat, D_profile_mat);
set(plot1,'LineStyle','none');
colorbar;
ylim([-5 5]);
xlim([0 4*pi]);
xlabel('\alpha\zeta');
ylabel('y');
title(sprintf('D1, Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));
subplot(5,3,14);
plot(D12_vec, xi_vec./sqrt(1-xi_vec.^2));
ylim([-5 5]);
ylabel('y');
title('D_0');
subplot(5,3,15)
hold on;
plot(abs(D_profile), xi_vec./sqrt(1-xi_vec.^2));
plot(real(D_profile), xi_vec./sqrt(1-xi_vec.^2), 'k--');
plot(imag(D_profile), xi_vec./sqrt(1-xi_vec.^2), 'k:');
hold off;
ylim([-5 5]);
ylabel('y');
ylabel('|D_1|');
title(sprintf('Re = %.4g, N = %d, alpha = %.4g, omega = %.4g + i%.4g', Rey, N, alpha, real(omega_vec(ii)), imag(omega_vec(ii))));