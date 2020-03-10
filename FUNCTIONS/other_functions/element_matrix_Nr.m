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
bh = dataElement.bh;

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
    dNU = diag(bh.dnu(Bmod));
%     [NU, dNU] = BHtoolNR (Bmod*1e3,2);
%     NU = diag(NU); dNU = diag(dNU);
else
    dNU = diag(Bmod*0);
end

D   = zeros([size(Dx) 2]); D(:,:,1) = Dx; D(:,:,2) = Dy;
Tnu = zeros([size(Dx) 4]);
            Tnu(:,:,1) = W*(NU + 2*dNU*byby); Tnu(:,:,3) = -W*2*dNU*bxby;
            Tnu(:,:,2) = -W*2*dNU*bybx;       Tnu(:,:,4) = W*(NU + 2*dNU*bxbx);
            
Tf  = zeros([size(Dx) 4]);
            Tf(:,:,1) = diag(Yeta(:));  Tf(:,:,3) = diag(-Xeta(:));
            Tf(:,:,2) = diag(-Ycsi(:)); Tf(:,:,4) = diag(Xcsi(:));
            
b = zeros([size(NU) 2]);
b(:,:,1) = Lxp; b(:,:,2) = Lyp;

DB = div(Tnu,D,0,1); DBF = div(Tf,DB,0,1);
Eij = -div(b,DBF,0,1);
end 