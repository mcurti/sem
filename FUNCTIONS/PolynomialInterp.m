%==========================================================================
% Current_Sources.m
% Created: 13.11.2017 - 13:11:32
% By: M. Curti
%
% Polynomial interpolation
%==========================================================================
function [ l ] = PolynomialInterp(x, xj)
%PolynomialInterp Polynomial interpolation matrix
%   x  = evaluation points
%   xj = node points
N = length(xj); M = length(x);
xj = reshape(xj,1,N);
wj = BarycentricWeights(xj);
l = zeros(M,N);

for ii = 1:M
    for jj = 1:N
        if x(ii) == xj(jj)
            l(ii,jj) = 1;
        else
            l(ii,jj) = wj(jj)/((x(ii)-xj(jj))*sum(wj./(x(ii)-xj)));
        end
    end
end

end

