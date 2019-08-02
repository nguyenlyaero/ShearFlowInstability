figure;
plot(alpha_vec, omega_i_max_vec);
ylim([0 inf]);
ylabel('\omega_i');
xlabel('\alpha');

figure;
plot(alpha_vec, omega_r_vec./alpha_vec);
ylabel('\omega_r/\alpha');
xlabel('\alpha');
