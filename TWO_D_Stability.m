%% Stability
N = 800;
alpha = 0.58;
[A,B] = TWO_D_CalcABTemporalCheb(alpha,N, paramArray);
[V,D] = eig(inv(B)*A);
omega_vec = diag(D);

% Sort
[~,I] = sort(imag(omega_vec));
omega_vec = omega_vec(I);
omega_vec = omega_vec(end:-1:1);
V = V(:,I); V = V(:,end:-1:1);

save('data_stability_800.mat', 'V', 'omega_vec', 'alpha');