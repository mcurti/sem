function [ out] = InterpolationPoints( csi, eta, csi_i, eta_i, s )
%InterpolationPoints Finds the interpolation points csi_i eta_i that
%corespond to Xm_i Ym_i nodes from FEM
%   Detailed explanation goes here

g1x = s{1,1}; g2x = s{1,2}; g3x = s{1,3}; g4x = s{1,4};
g1y = s{2,1}; g2y = s{2,2}; g3y = s{2,3}; g4y = s{2,4};

G1X = LagrangeInterp( csi_i, csi', g1x );
G1Y = LagrangeInterp( csi_i, csi', g1y );

G2X = LagrangeInterp( eta_i, eta', g2x );
G2Y = LagrangeInterp( eta_i, eta', g2y );

G3X = LagrangeInterp( csi_i, csi', g3x );
G3Y = LagrangeInterp( csi_i, csi', g3y );

G4X = LagrangeInterp( eta_i, eta', g4x );
G4Y = LagrangeInterp( eta_i, eta', g4y );
% Getting the mapped meshgrid

CSI = csi_i; ETA = eta_i;

Xm_i = 0.5*((1-CSI).*G4X + (1+CSI).*G2X + (1-ETA).*G1X + (1+ETA).*G3X) ...
    - 0.25*((1-CSI).*((1-ETA).*g1x(  1) + (1+ETA).*g3x(  1)) ...
          + (1+CSI).*((1-ETA).*g1x(end) + (1+ETA).*g3x(end)));

Ym_i = 0.5*((1-CSI).*G4Y + (1+CSI).*G2Y + (1-ETA).*G1Y + (1+ETA).*G3Y) ...
    - 0.25*((1-CSI).*((1-ETA).*g1y(  1) + (1+ETA).*g3y( 1)) ...
          + (1+CSI).*((1-ETA).*g1y(end) + (1+ETA).*g3y(end)));
      
      out = [Xm_i; Ym_i];
end

