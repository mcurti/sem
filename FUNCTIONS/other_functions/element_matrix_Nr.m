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
a    =  dataElement.a;
W = dataElement.W;
kk = dataElement.kk;


% A  = (Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
% B  = (Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
% C  = (Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));
W  = diag(W(:));


% Computing the Jacobian for each element
% 1. Computing the derivative matrices
Dx = diag(1./J(:))*( diag(Yeta(:))*Lx - diag(Ycsi(:))*Ly);
Dy = diag(1./J(:))*(-diag(Xeta(:))*Lx + diag(Xcsi(:))*Ly);
% F1 = A*Lx - B*Ly; F2 = C*Ly - B*Lx;

Bmod  = sqrt((Dx*a).^2 + (Dy*a).^2);
bx = diag(Dy*a); by = diag(-Dx*a);

bxbx = bx*bx; bxby = bx*by; bybx = by*bx; byby = by*by;

if max(abs(kk))>0
    bh = bh_class(2);
    dNU = diag(bh.dnu(Bmod));
%     [NU, dNU] = BHtoolNR (Bmod*1e3,2);
%     NU = diag(NU); dNU = diag(dNU);
else
    dNU = diag(Bmod*0);
end
% aTx = Lxp*W; aTy = Lyp*W;
% % Eij = -(aTx*NU*F1 + aTy*NU*F2 + 2*dNU*(aTx*bxbx + aTx*bxby)*F1 + (aTy*bybx + aTy*byby)*F2);
% % Eij = -(aTx*NU*F1 + aTy*NU*F2 + 2*(aTx*byby - aTx*bxby)*dNU*F1 + 2*(-aTy*bybx + aTy*bxbx)*dNU*F2);
% Eij = -(aTx*NU*F1 + aTy*NU*F2 + 2*aTx*(bxbx*dNU*F1 + bxby*dNU*F2) + 2*aTy*(bybx*dNU*F1 + byby*dNU*F2));

% a1 = ( diag(Yeta(:))*Lx - diag(Ycsi(:))*Ly);
% a2 = (-diag(Xeta(:))*Lx + diag(Xcsi(:))*Ly);
% 
% b1 = ( diag(Yeta(:))*Lxp - diag(Ycsi(:))*Lyp);
% b2 = (-diag(Xeta(:))*Lxp + diag(Xcsi(:))*Lyp);
% 
% JWdNU = diag(1./J(:))*W*2*dNU;
% Eij = -(b1*diag(1./J(:))*W*NU*a1 + ...
%         b2*diag(1./J(:))*W*NU*a2) + ...
%         b1*JWdNU*bxbx*a1    + b1*JWdNU*(bxby)*(a2) + ...
%         b2*JWdNU*(bybx)*a1 + b2*JWdNU*byby*(a2);
D   = zeros([size(Dx) 2]); D(:,:,1) = Dx; D(:,:,2) = Dy;
Tnu = zeros([size(Dx) 4]);
            Tnu(:,:,1) = W*(NU + 2*dNU*byby); Tnu(:,:,3) = -W*2*dNU*bxby;
            Tnu(:,:,2) = -W*2*dNU*bybx;       Tnu(:,:,4) = W*(NU + 2*dNU*bxbx);
            
Tf  = zeros([size(Dx) 4]);
            Tf(:,:,1) = diag(Yeta(:));  Tf(:,:,3) = diag(-Xeta(:));
            Tf(:,:,2) = diag(-Ycsi(:)); Tf(:,:,4) = diag(Xcsi(:));
            
b = zeros([size(NU) 2]);
b(:,:,1) = Lxp; b(:,:,2) = Lyp;

DB = div(Tnu,D); DBF = div(Tf,DB);
Eij = -div(b,DBF);
end 