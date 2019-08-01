function [h_ig_RT, cp_ig_R] = IdealEnthalpy(T, x, S1, S2)
% Function to return ideal gas enthalpy using NASA polynomials
    if nargin == 3
        coeffs1 = getNasaPoly(S1);
        h1_ig_RT = coeffs1(1) + coeffs1(2)*T/2 + coeffs1(3)*T^2/3 + coeffs1(4)*T^3/4 + coeffs1(5)*T^4/5 + coeffs1(6)/T;
        h_ig_RT = h1_ig_RT;
        
        cp_ig_R = coeffs1(1) + coeffs1(2)*T + coeffs1(3)*T^2 + coeffs1(4)*T^3 + coeffs1(5)*T^4;
        
    else % nargin ==4
        coeffs1 = getNasaPoly(S1);
        coeffs2 = getNasaPoly(S2);

        h1_ig_RT = coeffs1(1) + coeffs1(2)*T/2 + coeffs1(3)*T^2/3 + coeffs1(4)*T^3/4 + coeffs1(5)*T^4/5 + coeffs1(6)/T;
        h2_ig_RT = coeffs2(1) + coeffs2(2)*T/2 + coeffs2(3)*T^2/3 + coeffs2(4)*T^3/4 + coeffs2(5)*T^4/5 + coeffs2(6)/T;

        h_ig_RT = x*h1_ig_RT + (1-x)*h2_ig_RT;
        
        cp1_ig_R = coeffs1(1) + coeffs1(2)*T + coeffs1(3)*T^2 + coeffs1(4)*T^3 + coeffs1(5)*T^4;
        cp2_ig_R = coeffs2(1) + coeffs2(2)*T + coeffs2(3)*T^2 + coeffs2(4)*T^3 + coeffs2(5)*T^4;
        
        cp_ig_R = x*cp1_ig_R + (1-x)*cp2_ig_R;
    end
    
end