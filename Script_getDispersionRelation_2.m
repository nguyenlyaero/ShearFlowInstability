addpath("./DERIVESTsuite");
addpath("./PC-SAFT");
%p = parpool('local', 16);
%p.IdleTimeout = 120;
%% Baseflow
N = 800;
P_Pc = 1.18;
Ma = 0.5;
y_pb = 0.8;
paramArray = initParam(P_Pc, Ma, y_pb);
baseFlowArray = Calc_BaseFlow_Nodes(N, paramArray);
save Baseflow_800_0d8.mat baseFlowArray;
%% Stability
res = 100;
alpha_vec = linspace(0, 1.3, res)';
omega_i_max_vec = zeros(res,1);
omega_r_vec = zeros(res,1);
parfor i=1:res
    fprintf('%d \n', i);
    [omega_i_max_vec(i), omega_r_vec(i)] = getUnstableOmegaPos(alpha_vec(i), N, paramArray, baseFlowArray)
    fprintf('omega_i = %.4g \n omega_r = %.4g \n', omega_i_max_vec(i), omega_r_vec(i));
end

%figure;
%plot(alpha_vec, omega_i_max_vec);
%ylim([0 inf]);
%ylabel('\omega_i');
%xlabel('\alpha');
%title(sprintf('Re = %.4g, N = %d', Rey, N));

%figure;
%plot(alpha_vec, omega_r_vec./alpha_vec);
%ylabel('\omega_r/\alpha');
%xlabel('\alpha');
%title(sprintf('Re = %.4g, N = %d', Rey, N));

save('Dispersion_Pr_1d18_M_0d5_y_0d8_N_800.mat');
