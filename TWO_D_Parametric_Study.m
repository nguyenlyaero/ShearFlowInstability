addpath("../PC-SAFT")
Rey_vec = logspace(4, 5.6990, 5);
P_Pc_vec = linspace(1.1, 2,  10);
y_pb_vec = linspace(0, 1, 10);

%% Re1
for i = 3:3
    for j = 6:10
        if j == 1
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(1), y_pb_vec(j), 0.9, 1.0);
        else
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(1), y_pb_vec(j), alphaMax(i, j - 1), alphaMax(i, j - 1) - 0.05);
        end
        fprintf("##################################################### \n");
        fprintf("i = %d, j = %d, alphaMax = %.4g \n", i, j, alphaMax(i));
        fprintf("##################################################### \n");
        save("data1.mat");
    end
end

%% Re2
for i = 1:10
    for j = 1:10
        if j == 1
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(2), y_pb_vec(j), 0.9, 1.0);
        else
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(2), y_pb_vec(j), alphaMax(i, j - 1), alphaMax(i, j - 1) - 0.05);
        end
        fprintf("##################################################### \n");
        fprintf("i = %d, j = %d, alphaMax = %.4g \n", i, j, alphaMax(i));
        fprintf("##################################################### \n");
        save("data2.mat");
    end
end

%% Re3
for i = 1:1
    for j = 9:10
        if j == 1
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(3), y_pb_vec(j), 0.9, 1.0);
        else
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(3), y_pb_vec(j), alphaMax(i, j - 1), alphaMax(i, j - 1) - 0.05);
        end
        fprintf("##################################################### \n");
        fprintf("i = %d, j = %d, alphaMax = %.4g \n", i, j, alphaMax(i));
        fprintf("##################################################### \n");
        save("data3.mat");
    end
end

%% Re4
for i = 1:10
    for j = 1:10
        if i == 1
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(4), y_pb_vec(j), 0.9, 1.0);
        else
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(4), y_pb_vec(j), alphaMax(i, j - 1), alphaMax(i, j - 1) - 0.05);
        end
        fprintf("##################################################### \n");
        fprintf("i = %d, j = %d, alphaMax = %.4g \n", i, j, alphaMax(i));
        fprintf("##################################################### \n");
        save("data4.mat");
    end
end

%% Re5
for i = 1:10
    for j = 1:10
        if i == 1
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(5), y_pb_vec(j), 0.9, 1.0);
        else
            alphaMax(i,j) = parametricStudy(P_Pc_vec(i), Rey_vec(5), y_pb_vec(j), alphaMax(i, j - 1), alphaMax(i, j - 1) - 0.05);
        end
        fprintf("##################################################### \n");
        fprintf("i = %d, j = %d, alphaMax = %.4g \n", i, j, alphaMax(i));
        fprintf("##################################################### \n");
        save("data5.mat");
    end
end