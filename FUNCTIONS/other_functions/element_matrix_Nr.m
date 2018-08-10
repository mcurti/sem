%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function Eij = element_matrix_Nr(dataElement)
Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta = dataElement.Xeta;  Yeta = dataElement.Yeta;
J    = dataElement.J;
NU   = diag(dataElement.nu(:));
a   =  dataElement.a;
W = dataElement.W;
kk = dataElement.kk;


A  = (Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B  = (Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C  = (Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));
W  = diag(W(:));


% Computing the Jacobian for each element
% 1. Computing the derivative matrices
Dx = diag(1./J(:))*( diag(Yeta(:))*Lx - diag(Ycsi(:))*Ly);
Dy = diag(1./J(:))*(-diag(Xeta(:))*Lx + diag(Xcsi(:))*Ly);
F1 = A*Lx - B*Ly; F2 = C*Ly - B*Lx;

Bmod  = sqrt((Dx*a).^2 + (Dy*a).^2);
by = diag(Dy*a); bx = diag(Dx*a);

bxbx = bx.*bx; bxby = bx.*by; bybx = by.*bx; byby = by.*by;

if max(abs(kk))>0
    bh = bh_class(2);
    dNU = diag(bh.dnu(Bmod*1e3));
else
    dNU = diag(Bmod*0);
end

Eij = -(Lxp*W*((NU + 2*dNU*byby)*F1     +      -2*dNU*bxby*F2) ...
       +Lyp*W*(       -2*dNU*bybx*F1    +   (NU + 2*dNU*bxbx)*F2));

end 