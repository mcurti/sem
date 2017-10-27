function [ Xcsi, Ycsi, Xeta, Yeta, J ] = TransfiniteQuadMetrics1( csi, eta, s )
%TransfiniteQuadMap Returns the mapped meshgrid according to the boundary
%curves
%   csi      - x points
%   eta      - y points
%   curves_x - points of the 4 curves for x output
%   curves_y - points of the 4 curves for y output

% reordering the input data
csi = reshape(csi,length(csi),1); eta = reshape(eta,length(eta),1);
% curves_x = reshape(curves_x,4,length(csi));
% curves_y = reshape(curves_y,4,length(eta));

[CSI, ETA] = meshgrid(csi, eta);
% Getting the mapped curves:

G1x = s{1,1}; G2x = s{1,2}; G3x = s{1,3}; G4x = s{1,4};
G1y = s{2,1}; G2y = s{2,2}; G3y = s{2,3}; G4y = s{2,4};

[G1X, G2Y] = meshgrid(G1x, G2y); [G3X, G4Y] = meshgrid(G3x, G4y);
[G1Y, G2X] = meshgrid(G1y, G2x); [G3Y, G4X] = meshgrid(G3y, G4x);

% [G1x, G1y, G1X, G1Y]   = EvaluateAt(csi,curves_x(1,:),curves_y(1,:),'csi');
% [~  , ~  , G2X, G2Y]   = EvaluateAt(eta,curves_x(2,:),curves_y(2,:),'eta');
% [G3x, G3y, G3X, G3Y]   = EvaluateAt(csi,curves_x(3,:),curves_y(3,:),'csi');
% [~  , ~  , G4X, G4Y]   = EvaluateAt(eta,curves_x(4,:),curves_y(4,:),'eta');

G1px     = LagrangeInterpDerivative(csi,csi',G1x);
G2px     = LagrangeInterpDerivative(eta,eta',G2x);
G3px     = LagrangeInterpDerivative(csi,csi',G3x);
G4px     = LagrangeInterpDerivative(eta,eta',G4x);

G1py     = LagrangeInterpDerivative(csi,csi',G1y);
G2py     = LagrangeInterpDerivative(eta,eta',G2y);
G3py     = LagrangeInterpDerivative(csi,csi',G3y);
G4py     = LagrangeInterpDerivative(eta,eta',G4y);

[G1pX, G2pY] = meshgrid(G1px, G2py); [G3pX, G4pY] = meshgrid(G3px, G4py);
[G1pY, G2pX] = meshgrid(G1py, G2px); [G3pY, G4pX] = meshgrid(G3py, G4px);

% [~, ~, G2pX, G2pY]     = DerivativeAt(eta,curves_x(2,:),curves_y(2,:),'eta');
% [~, ~, G3pX, G3pY]     = DerivativeAt(csi,curves_x(3,:),curves_y(3,:),'csi');
% [~, ~, G4pX, G4pY]     = DerivativeAt(eta,curves_x(4,:),curves_y(4,:),'eta');
% Getting the metrics and Jacobian
% G4pY = -G4pY; 
Xcsi = 0.5*(G2X - G4X + (1-ETA).*G1pX + (1+ETA).*G3pX) - ...
       0.25*((1-ETA)*(G1x(end)-G1x(1)) + (1+ETA)*(G3x(end)-G3x(1)));
Ycsi = 0.5*(G2Y - G4Y + (1-ETA).*G1pY + (1+ETA).*G3pY) - ...
       0.25*((1-ETA)*(G1y(end)-G1y(1)) + (1+ETA)*(G3y(end)-G3y(1)));

Xeta = 0.5*((1-CSI).*G4pX + (1+CSI).*G2pX + G3X - G1X) - ...
       0.25*((1-CSI)*(G3x(1)-G1x(1)) + (1+CSI)*(G3x(end)-G1x(end)));
Yeta = 0.5*((1-CSI).*G4pY + (1+CSI).*G2pY + G3Y - G1Y) - ...
       0.25*((1-CSI)*(G3y(1)-G1y(1)) + (1+CSI)*(G3y(end)-G1y(end)));
J    = Xcsi.*Yeta - Ycsi.*Xeta;

% J = J'; Xcsi = Xcsi'; Ycsi = Ycsi'; Xeta = Xeta';
% Yeta = Yeta';
end

