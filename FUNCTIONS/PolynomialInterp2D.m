function [ lxy ] = PolynomialInterp2D( x, y, Xj, Yj)
%PolynomialInterp2D returns the matrix for 2D interpolation
%   Detailed explanation goes here
%   Works only for rectangular domains

[M, N]  = size(Xj); L = numel(x);
wx = zeros(M,N); wy = zeros(M,N); tmp = zeros(M,N); lxy = zeros(L,N*M);

for k = 1:M
    wx(k,:) = BarycentricWeights(Xj(k,:));
end
% wx
for k = 1:N
    wy(:,k) = BarycentricWeights(Yj(:,k));
end
% wy
% wx

%     assignin('base', 'wy', wy);
%     assignin('base', 'wx', wx);
%     assignin('base', 'Yj', Yj);
%     assignin('base', 'Xj', Xj);
for k = 1:L
    for ii = 1:M
        for jj = 1:N
            if x(k) == Xj(ii,jj) && y(k) == Yj(ii,jj)
               tmp(ii,jj) = 1;
            elseif x(k) == Xj(ii,jj)
               tmp(ii,jj) = 1 * ...
                            wy(ii,jj)./((y(k)-Yj(ii,jj)).*sum(wy(:,jj)./(y(k)-Yj(:,jj))));
            elseif y(k) == Yj(ii,jj)
               tmp(ii,jj) = wx(ii,jj)./((x(k)-Xj(ii,jj)).*sum(wx(ii,:)./(x(k)-Xj(ii,:)))) * ...
                            1;
            else
               tmp(ii,jj) = wx(ii,jj)./((x(k)-Xj(ii,jj)).*sum(wx(ii,:)./(x(k)-Xj(ii,:)))) * ...
                            wy(ii,jj)./((y(k)-Yj(ii,jj)).*sum(wy(:,jj)./(y(k)-Yj(:,jj))));
            end
        end
    end
    lxy(k,:) = tmp(:);
end
end

