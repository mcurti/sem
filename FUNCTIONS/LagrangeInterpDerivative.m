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
% selecting non equal points
% xii = true(size(xj));

xj = reshape(xj,length(xj),1);
fj = reshape(fj,length(fj),1);

wj = reshape(wj,size(xj));

Px  = LagrangeInterp(x, xj, fj);

for k = 1:length(x)
    if sum(x(k)==xj) == 0
        % Not at the nodes
        Ppx(k) = sum(wj./(x(k) - xj).*(Px(k)-fj)./(x(k) - xj))./sum(wj./(x(k) - xj));
    else
        % At the nodes
        xjj = not(x(k)==xj); xii = not(xjj);
        Ppx(k) = -1./wj(xii).*sum(wj(xjj).*(Px(k)-fj(xjj))./(x(k) - xj(xjj)));
    end
end

end

