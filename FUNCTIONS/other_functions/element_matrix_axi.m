%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function Eij = element_matrix_axi(dataElement)
Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta = dataElement.Xeta;  Yeta = dataElement.Yeta;
Xminv   = diag(1./dataElement.Xm(:));
J    = dataElement.J;
nu   = dataElement.nu;
W = dataElement.W;


A = (Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B = (Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C = (Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));


NU = Xminv*diag(nu(:));
W = diag(W(:));

F = zeros([size(NU) 2]);
b = F;
F(:,:,1) = W*NU*(A*Lx - B*Ly); F(:,:,2) = W*NU*(C*Ly - B*Lx);
b(:,:,1) = Lxp; b(:,:,2) = Lyp;

Eij = - div(b,F,1,1);
end 