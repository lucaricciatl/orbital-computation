

function Tzxz = Tzxz(alpha,beta,gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tz1 =  T_matrix.Tz(alpha);
Tx2 = T_matrix.Tx(beta);
Tz3 = T_matrix.Tz(gamma);
Tzxz = Tz1*Tx2*Tz3;
end

