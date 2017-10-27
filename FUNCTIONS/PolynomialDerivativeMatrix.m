function [ D ] = PolynomialDerivativeMatrix( x )
%PolynomoalDerivativeMatrix First Derivative Approximation Matrix algorithm
%37
%   Compute the first derivative approximation matrix in points "x_j"

N = length(x);

% Computing the barycentric weights
w = BarycentricWeights(x);

% Computing the derivative matrix

D = zeros(N);

for ii = 1:N
    D(ii,ii) = 0;
    for jj = 1:N
        if jj ~= ii
            D(ii,jj) = w(jj)/(w(ii)*(x(ii)-x(jj)));
            D(ii,ii) = D(ii,ii) - D(ii,jj);
        end
    end
end
end

