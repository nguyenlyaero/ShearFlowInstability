function [k_bin] = getBinaryInteraction(S1, S2)
%GETBINARYINTERACTION Summary of this function goes here
%   Detailed explanation goes here
    
    % Delters 1993
    k_H2_O2 = 1-0.92683;
    
    % Lease Squares with CoolProp
    k_CH4_O2 = 0.098561433743562;
    
    %
    k_CH4_Ethane = 0;
    k_CH4_N2 = 0;
    
    switch 1
        case strcmp(S1,'H2') && strcmp(S2,'O2')
            k_bin = k_H2_O2;
        case strcmp(S1,'CH4') && strcmp(S2,'O2')
            k_bin = k_CH4_O2;
        case strcmp(S1,'CH4') && strcmp(S2,'N2')
            k_bin = k_CH4_N2;
        case strcmp(S1,'CH4') && strcmp(S2,'C2H4')
            k_bin = k_CH4_Ethane;
        otherwise
            k_bin = 0;
    end
end

