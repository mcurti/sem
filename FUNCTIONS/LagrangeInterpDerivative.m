function [ Ppx ] = LagrangeInterpDerivative( x, xj, fj )
%LAGRANGEINTERPDERIVATIVE  Performes Lagrange interpolation derivative
%
%   x  = list of points where the interpolation is evaluated
%   xj = list of points where the fj is given

% computing the barycentric weights
% N = length(xj);

wj = BarycentricWeights(xj);
% evaluating the function fj on x nodes
Ppx = zeros(size(x));
Px  = LagrangeInterp(x, xj, fj);
% selecting non equal points
% xii = true(size(xj));


for k = 1:length(x)
    if sum(x(k)==xj) == 0
    Ppx(k) = sum(wj.*(Px(k)-fj)./(x(k) - xj).^2)./sum(wj./(x(k) - xj));
    else
        xii = not(x(k)==xj);
        Ppx(k) = -1./wj(k).*sum(wj(xii).*(fj(k)-fj(xii))./(x(k) - xj(xii)));
    end
end

end

