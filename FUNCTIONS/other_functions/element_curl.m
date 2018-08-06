%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function Eij = element_curl(dataElement)
Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta = dataElement.Xeta;  Yeta = dataElement.Yeta;
J    = dataElement.J;     PHI = dataElement.PHI;
% nu   = dataElement.nu;
% W = dataElement.W;


A = (Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B = (Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C = (Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));

% W = diag(W(:));


Eij =   reshape(Lx*(A*(Lx*PHI(:)) - B*(Ly*PHI(:))) + Ly*(C*(Ly*PHI(:)) - B*(Lx*PHI(:))),size(J))./J;

end 