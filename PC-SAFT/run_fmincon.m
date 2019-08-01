function k_CH4_O2 = run_fmincon()
    import py.CoolProp.CoolProp.PropsSI;
    P = 1.18*50.43*1e5; % Pa
    T = 201.4; % K
    x_vec = linspace(0,1,100)';
    rho_vec_CoolProp = zeros(100,1);
    for i=1:100
        string = sprintf('CH4[%.8g]&O2[%.8g]', x_vec(i), 1-x_vec(i));
        rho_vec_CoolProp(i) = PropsSI('D', 'T', T, 'P', P, string);
    end

    k_CH4_O2 = fmincon(@Objective_Fun,0,[1;-1],[1;0]);

    function out = Objective_Fun(k_CH4_O2)
        rho_vec_PC_SAFT = zeros(100,1);
        for i=1:100
            rho_hat = PC_SAFT_PT(P,T,x_vec(i),'CH4&O2',k_CH4_O2);
            W = 16*x_vec(i) + 32*(1-x_vec(i));
            rho_vec_PC_SAFT(i) = rho_hat*W;
        end
        out = sum(((rho_vec_PC_SAFT-rho_vec_CoolProp)./rho_vec_CoolProp).^2);
        % figure;
        % hold on;
        % plot(x_vec, rho_vec_PC_SAFT);
        % plot(x_vec, rho_vec_CoolProp);
        % hold off;
        % 
        % figure;
        % plot(x_vec, (rho_vec_PC_SAFT-rho_vec_CoolProp)./rho_vec_CoolProp);
    end
end

