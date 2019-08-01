function [Tc, Pc, acentric, M, m , epsilon_k, sigma, sigma_coll, dipole, Vc] = getSpecies(S)
    % Species data
    
    switch 1
        case strcmp(S, 'H2')
            % Senol 2011
            Tc = 33.2; % K
            Pc = 13.0e5; % Pa
            acentric = -0.218;
            M = 2; % [g/mol]
            m = 0.935864;
            epsilon_k = 25.62934; % [K]
            sigma = 2.912599; % [Angstrom]
        case strcmp(S, 'O2')
            % Stoll Vrabec Hasse 2003
            Tc = 154.6; % K
            Pc = 50.4e5; % Pa
            acentric = 0.025;
            M = 32.00; % [g/mol]
            m = 1.1457;
            epsilon_k = 113.98; % [K]
            sigma = 3.1711; % [Angstrom]
            sigma_coll = 3.54; % [Angstrom] (Bird)
            dipole = 0; %  [D] Dipole momment
            Vc = 0.0734; % [l/mol] Critical molar specific volume
        case strcmp(S, 'N2')
            % Stoll Vrabec Hasse 2003
            Tc = 126.2; % K
            Pc = 33.9e5; % Pa
            acentric = 0.039;
            M = 28; % [g/mol]
            m = 1.2365;
            epsilon_k = 89.492; % [K]
            sigma = 3.2975; % [Angstrom]
            sigma_coll = 3.54; % [Angstrom] PLACEHOLDER
            dipole = 0;
            Vc = 0.0897; % [l/mol] Critical molar specific volume
        case strcmp(S, 'CH4')
            % Gross Sadowski 2001
            Tc = 190.6; % K
            Pc = 46.0e5; % Pa
            acentric = 0.008;
            M = 16.043; % [g/mol]
            m = 1; % # Segments per chain
            sigma = 3.7039; % Segment Diameter [Angstrom]
            epsilon_k = 150.03; % Potential Depth [K]
            sigma_coll = 3.80; %[Angstrom] Molecular collision diamater (Bird)
            dipole = 0; %  [D] Dipole momment
            Vc = 0.0991; % [l/mol] Critical molar specific volume
        case strcmp(S, 'CO2')
            % Gross Sadowski 2001
            Tc = 304.2; % K
            Pc = 73.8e5; % Pa
            acentric = 0.224;
            M = 44.01; % [g/mol]
            m = 2.0729; % # Segments per chain
            sigma = 2.7852; % Segment Diameter [Angstrom]
            epsilon_k=169.21; % Potential Depth [K]
            sigma_coll = 3.54; % [Angstrom] PLACEHOLDER
            dipole = 0;
            Vc = 0.09412; % [l/mol] Critical molar specific volume
        case strcmp(S, 'C2H4')
            % Gross Sadowski 2001
            Tc = 305.4; % K
            Pc = 48.8e5; % Pa
            acentric = 0.099;
            M = 30.07; % [g/mol]
            m = 1.6069; % # Segments per chain
            sigma = 3.5206; % Segment Diameter [Angstrom]
            epsilon_k = 191.42; % Potential Depth [K]
        case strcmp(S, 'C4H10')
            % Gross Sadowski 2001
            Tc = 425.2;
            Pc = 38.0e5;
            acentric =  0.199;
            M = 58.123;           
            m = 2.3316;
            sigma = 3.7086;
            epsilon_k = 222.88;
    end
end

