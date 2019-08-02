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

    omega_best = omega_vec(1);    
    omega_i = imag(omega_vec(1));
    omega_r = real(omega_vec(1));
    
    if omega_i > 0
        lowestPhaseSpeed = 0.8;
        for i=1:length(omega_vec)
            if imag(omega_vec(i)) < 0
                break;
            elseif abs(real(omega_vec(i))/alpha) < lowestPhaseSpeed
                omega_best = omega_vec(i);
                lowestPhaseSpeed = abs(real(omega_vec(i))/alpha);
            end
        end
        if lowestPhaseSpeed == 0.8
            omega_best = omega_vec(i);    
        end
        omega_i = imag(omega_best);
        omega_r = real(omega_best);
    end
end

