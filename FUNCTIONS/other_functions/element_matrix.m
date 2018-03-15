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


A = nu.*(Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B = nu.*(Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C = nu.*(Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));

W = diag(W(:));


Eij =   - (Lxp*W*(A*Lx - B*Ly) + Lyp*W*(C*Ly - B*Lx));

end 