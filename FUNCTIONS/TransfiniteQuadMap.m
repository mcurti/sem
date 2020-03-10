function [ X, Y ] = TransfiniteQuadMap( csi, eta, s )
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
% [G1x, G1y, G1X, G1Y]     = EvaluateAt1(s{1,1},s{2,1},'csi');
% [~, ~, G2X, G2Y]     = EvaluateAt1(s{1,2},s{2,2},'eta');
% [G3x, G3y, G3X, G3Y]     = EvaluateAt1(s{1,3},s{2,3},'csi');
% [~, ~, G4X, G4Y]     = EvaluateAt1(s{1,4},s{2,4},'eta');

% Getting the mapped meshgrid

X = 0.5*((1-CSI).*G4X + (1+CSI).*G2X + (1-ETA).*G1X + (1+ETA).*G3X) ...
    - 0.25*((1-CSI).*((1-ETA).*G1x(  1) + (1+ETA).*G3x(  1)) ...
          + (1+CSI).*((1-ETA).*G1x(end) + (1+ETA).*G3x(end)));

Y = 0.5*((1-CSI).*G4Y + (1+CSI).*G2Y + (1-ETA).*G1Y + (1+ETA).*G3Y) ...
    - 0.25*((1-CSI).*((1-ETA).*G1y(  1) + (1+ETA).*G3y( 1)) ...
          + (1+CSI).*((1-ETA).*G1y(end) + (1+ETA).*G3y(end)));

end

