%==========================================================================
% hbnp_r.m
% Created: 09.02.2018 - 11:50:36
% By: M. Curti
%
% the set of linearly independent functions for polar coordinate system
%==========================================================================
function [h_p, h_n, b_p, b_n] = hbnp_r(w_n,r,s1,s2)
rk1 = r/s1; rk2 = r/s2;

% h_p = r.^(w_n-1)./(s1.^w_n).*w_n;
% h_n = r.^(-w_n-1)./(s2.^(-w_n)).*(-w_n);
h_p = rk1.^(w_n-1).*w_n/s1; h_n = rk2.^(-w_n-1).*(-w_n)/s2;

b_p = rk1.^(w_n);   b_n = rk2.^(-w_n);
end
