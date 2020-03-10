%==========================================================================
% np_exp.m
% Created: 09.02.2018 - 11:50:36
% By: M. Curti
%
% the set of linearly independent functions
%==========================================================================
function [p_exp, n_exp] = np_exp(w_n,y,s1,s2)
    p_exp = exp( w_n*(y-s1));
    n_exp = exp(-w_n*(y-s2));
end

