function [ w ] = BarycentricWeights( x )
%BarycentricWeights Is returning the weights for Lagrange interpolation
%algorithm 30
%   For a set of points "x_i" on x axis the functions returns the weights
%   for Lagrange interpolation

N = length(x);

w = ones(1,N);

for jj = 2:(N)
    for k = 1:(jj-1)
           w(k)  = w(k) *(x(k ) - x(jj));
           w(jj) = w(jj)*(x(jj) - x(k ));
    end
end

w = 1./w;