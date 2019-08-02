addpath("./DERIVESTsuite");
addpath("./PC-SAFT");
%% Stability
N = 800;
P_Pc = 1.18;
Ma = 0.5;
y_pb = 0.8;
paramArray = initParam(P_Pc, Ma, y_pb);
% baseFlowArray = Calc_BaseFlow_Nodes(N, paramArray);
% save Baseflow_2000_0.mat baseFlowArray;
load Baseflow_800_0d8.mat;
%%
alpha = 0.52;
[A,B] = TWO_D_CalcABTemporalCheb(alpha,N, paramArray, baseFlowArray);
[V,D] = eig(inv(B)*A);
omega_vec = diag(D);

% Sort
[~,I] = sort(imag(omega_vec));
omega_vec = omega_vec(I);
omega_vec = omega_vec(end:-1:1);
V = V(:,I); V = V(:,end:-1:1);

save('data_stability__Pr_1d18_M_0d5_y_0d8_N_800_alpha_0d52.mat', 'V', 'omega_vec', 'alpha');