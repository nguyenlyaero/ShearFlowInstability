function [s_ig_R] = IdealEntropy(T, P, x, S1, S2)
% Function to return ideal gas entropy using NASA polynomials
% T: [K]
% P: [Pa]
    P0 = 1e5; %[Pa] Reference pressure 

    if nargin == 4
        coeffs1 = getNasaPoly(S1);
        
        s_ig_R = coeffs1(1)*log(T) + coeffs1(2)*T + coeffs1(3)*T^2/2 + coeffs1(4)*T^3/3 + coeffs1(5)*T^4/4 + coeffs1(7);
        s_ig_R = s_ig_R - log(P/P0);
        
    else % nargin == 5
        coeffs1 = getNasaPoly(S1);
        coeffs2 = getNasaPoly(S2);
        
        s1_ig_R = coeffs1(1)*log(T) + coeffs1(2)*T + coeffs1(3)*T^2/2 + coeffs1(4)*T^3/3 + coeffs1(5)*T^4/4 + coeffs1(7);
        s2_ig_R = coeffs2(1)*log(T) + coeffs2(2)*T + coeffs2(3)*T^2/2 + coeffs2(4)*T^3/3 + coeffs2(5)*T^4/4 + coeffs2(7);
        
        s_ig_R = x*s1_ig_R + (1-x)*s2_ig_R - x*log(x) - (1-x)*log(1-x);
        s_ig_R = s_ig_R - log(P/P0);
        
    end
    
end