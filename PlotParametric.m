load data1.mat;

figure;
hold on;
plot(y_pb_vec, alphaMax(1,:), 'o-');
plot(y_pb_vec, alphaMax(2,:), 'o-');
plot(y_pb_vec, alphaMax(4,:), 'o-');
plot(y_pb_vec, alphaMax(5,:), 'o-');
plot(y_pb_vec, alphaMax(7,:), 'o-');
plot(y_pb_vec, alphaMax(8,:), 'o-');
plot(y_pb_vec, alphaMax(9,:), 'o-');
hold off;
xlabel('y_{pb}');
ylabel('\alpha_{max}');
xlim([0 0.9]);
title('Re = 1E4');


figure;
hold on;
plot(y_pb_vec, alphaMax(1,:), 'o-');
plot(y_pb_vec([1:2 4:8 10]), alphaMax(2,[1:2 4:8 10]), 'o-');
plot(y_pb_vec(1:9), alphaMax(4,1:9), 'o-');
plot(y_pb_vec([1:2 4:5 7 9]), alphaMax(5,[1:2 4:5 7 9]), 'o-');
plot(y_pb_vec([1:4 6 8 9]), alphaMax(7,[1:4 6 8 9]), 'o-');
hold off;
legend('P_r = 1.1', 'P_r = 1.2', 'P_r = 1.4', 'P_r = 1.5', 'P_r = 1.7');
xlabel('y_{pb}');
ylabel('\alpha_{max}');
xlim([0 0.9]);
title('Re = 1E4');

load data2.mat;

figure;
hold on;
plot(y_pb_vec, alphaMax(1,:), 'o-');
plot(y_pb_vec([1:2 4:8 10]), alphaMax(2,[1:2 4:8 10]), 'o-');
plot(y_pb_vec(1:9), alphaMax(4,1:9), 'o-');
plot(y_pb_vec([1:2 4:5 7 9]), alphaMax(5,[1:2 4:5 7 9]), 'o-');
plot(y_pb_vec([1:4 6 8 9]), alphaMax(7,[1:4 6 8 9]), 'o-');
hold off;
legend('P_r = 1.1', 'P_r = 1.2', 'P_r = 1.4', 'P_r = 1.5', 'P_r = 1.7');
xlabel('y_{pb}');
ylabel('\alpha_{max}');
xlim([0 0.9]);
title('Re = 2.6E4');

load data3.mat;

figure;
hold on;
plot(y_pb_vec, alphaMax(1,:), 'o-');
plot(y_pb_vec([1:2 4:8 10]), alphaMax(2,[1:2 4:8 10]), 'o-');
plot(y_pb_vec(1:9), alphaMax(4,1:9), 'o-');
plot(y_pb_vec([1:2 4:5 7 9]), alphaMax(5,[1:2 4:5 7 9]), 'o-');
plot(y_pb_vec([1:4 6 8 9]), alphaMax(7,[1:4 6 8 9]), 'o-');
hold off;
legend('P_r = 1.1', 'P_r = 1.2', 'P_r = 1.4', 'P_r = 1.5', 'P_r = 1.7');
xlabel('y_{pb}');
ylabel('\alpha_{max}');
xlim([0 0.9]);
title('Re = 7E4');