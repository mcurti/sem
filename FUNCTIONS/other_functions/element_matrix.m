%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function Eij = element_matrix(dataElement)
Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta = dataElement.Xeta;  Yeta = dataElement.Yeta;
J    = dataElement.J;
nu   = dataElement.nu;
W = dataElement.W;


A = (Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B = (Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C = (Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));


NU = diag(nu(:));
W = diag(W(:));

F = zeros([size(NU) 2]);
b = F;
F(:,:,1) = W*NU*(A*Lx - B*Ly); F(:,:,2) = W*NU*(C*Ly - B*Lx);
b(:,:,1) = Lxp; b(:,:,2) = Lyp;

% a1 = ( diag(Yeta(:))*Lx - diag(Ycsi(:))*Ly);
% a2 = (-diag(Xeta(:))*Lx + diag(Xcsi(:))*Ly);
% 
% b1 = ( diag(Yeta(:))*Lxp - diag(Ycsi(:))*Lyp);
% b2 = (-diag(Xeta(:))*Lxp + diag(Xcsi(:))*Lyp);
% 
% Eij = -(b1*NU*diag(1./J(:))*W*a1 + b2*NU*diag(1./J(:))*W*a2);
% Eij =   - (Lxp*W*NU*(A*Lx - B*Ly) + Lyp*W*NU*(C*Ly - B*Lx));

Eij = - div(b,F);
end 