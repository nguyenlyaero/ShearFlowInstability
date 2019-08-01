function omega = getUnstableOmega(alpha, N, paramArray, baseFlowArray)
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
    
    omega = imag(omega_vec(1));
end
