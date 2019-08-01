N = 800;
[xi_vec, ~, ~, ~, ~] = Dmat(N);
load(['Base_Flow_' sprintf('%d', N) '.mat']);

cp_vec = zeros(N,1);
for i=1:N
    x_CH4 = M_O2*Y0_vec(i)/(M_CH4*(1-Y0_vec(i))+M_O2*Y0_vec(i));
    M = x_CH4*M_CH4 + (1-x_CH4)*M_O2;
    R = Rbar/M*1000;
    rho_hat = PC_SAFT_PT_Cubic(P_ref, T0_vec(i)*T_ref, x_CH4, 'CH4', 'O2', 0.5);
    cp_vec(i) = PC_SAFT_DT_cp(rho_hat, T0_vec(i)*T_ref, x_CH4, 'CH4', 'O2')*R/cp_ref;
end

figure;
hold on;
plot(xi_vec, u0__x_vec);
plot(xi_vec, P_vec);
plot(xi_vec, T0_vec*T_ref/T_base);
plot(xi_vec, Y0_vec);
hold off;
xlabel('\xi');
legend('u_0/U_c', 'P_0/P_c', 'T_0/T_c', 'Y_{CH4}');
title('P_{ctr}/P_{c,O2}=2.00');


figure;
hold on;
plot(xi_vec, mu_vec);
plot(xi_vec, lambda_vec);
plot(xi_vec, D12_vec);
plot(xi_vec, cp_vec);
hold off;
xlabel('\xi');
legend('\mu_0/\mu_c', '\lambda_0/lambda_c', 'D_0/D_c', 'cp_0/cp_{ref}');
title('P_{ctr}/P_{c,O2}=2.00');