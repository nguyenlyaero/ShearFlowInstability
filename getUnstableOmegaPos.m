function [omega_i, omega_r] = getUnstableOmegaPos(alpha, N, paramArray, baseFlowArray)
    [A,B] = TWO_D_CalcABTemporalCheb(alpha,N, paramArray, baseFlowArray);
    [~,D] = eig(B\A);

    omega_vec = diag(D);

    % Filter
    if alpha~=0
        I = find(abs(real(omega_vec) - alpha)<0.5*alpha);
    else
        I = find(abs(real(omega_vec) - alpha)<0.05);
    end
    omega_vec = omega_vec(I);
    
    % Sort
    [~,I] = sort(imag(omega_vec));
    omega_vec = omega_vec(I);
    omega_vec = omega_vec(end:-1:1);
    
    omega_i = imag(omega_vec(1));
    omega_r = real(omega_vec(1));
    
    if omega_i > 0
        diff = 100;
        for i=1:length(omega_vec)
            if imag(omega_vec(i)) > 0 && abs(real(omega_vec(i)) - alpha) < diff
                omega_best = omega_vec(i);
                diff = abs(real(omega_vec(i)) - alpha);
            end
        end
        omega_i = imag(omega_best);
        omega_r = real(omega_best);
    end
end

