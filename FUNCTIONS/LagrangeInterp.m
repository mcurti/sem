function [ Px ] = LagrangeInterp( x, xj, fj )
%LAGRANGE INTERPOLATION  Performes Lagrange interpolation on points x
%   x  = list of points where the interpolation is evaluated
%   xj = list of points where the fj is given

xj = reshape(xj,length(xj),1);
fj = reshape(fj,length(fj),1);
% computing the barycentric weights
wj = BarycentricWeights(xj);

wj = reshape(wj,size(xj));

% evaluating the function fj on x nodes
Px = zeros(size(x));

for k = 1:length(x)
    if sum(x(k)==xj) == 0
    Px(k) = sum(fj.*wj./(x(k) - xj))./sum(wj./(x(k) - xj));
    else
        a = fj(x(k)==xj);
        Px(k) = a(1);
    end
end

end

