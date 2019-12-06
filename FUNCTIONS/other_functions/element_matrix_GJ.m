%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function [Eij, Jij] = element_matrix_GJ(dataElement)
% Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Dx  = dataElement.Dx;
Dy  = dataElement.Dy;
Xcsi = diag(dataElement.Xcsi(:));  Ycsi = diag(dataElement.Ycsi(:));
Xeta = diag(dataElement.Xeta(:));  Yeta = diag(dataElement.Yeta(:));
% Jin    = diag(1./dataElement.J(:));
NU   = dataElement.nu(:);
a    =  dataElement.a;
W    = dataElement.W(:);
kk = dataElement.kk;
bh = dataElement.bh;


%% Computing the magnetic Field and reluctances
% 1. Computing the derivative matrices




bx = Dy*a; by = -Dx*a;

Bmod  = sqrt(by.^2 + bx.^2);

bxbx = 2*bx.*bx; bxby = 2*bx.*by; byby = 2*by.*by;

if max(abs(kk))>0
    dNU = bh.dnu(Bmod);
else
    dNU = zeros(numel(Bmod),1);
end

%% Creating the matrices
% Derivative matrices
D   = zeros([size(Dx) 2]); D(:,:,1) = Dx; D(:,:,2) = Dy;

% Reluctance tensors
WNU = W.*NU; WdNU = W.*dNU;
Tnu = zeros([size(Dx) 4]);
         Tnu(:,:,1) = diag(WNU + WdNU.*byby); Tnu(:,:,3) = diag(-WdNU.*bxby);
         Tnu(:,:,2) = Tnu(:,:,3); Tnu(:,:,4) = diag(WNU + WdNU.*bxbx);
Tnu0 = zeros([size(Dx) 4]);
            Tnu0(:,:,1) = diag(WNU);
                                         Tnu0(:,:,4) = diag(WNU);

% Flux tensors
Tf  = zeros([size(Dx) 4]);
            Tf(:,:,1) = Yeta;  Tf(:,:,3) = -Xeta;
            Tf(:,:,2) = -Ycsi; Tf(:,:,4) = Xcsi;
            
b = zeros([size(Dx) 2]);
b(:,:,1) = Lxp; b(:,:,2) = Lyp;

DB = div(Tnu,D,0,1); DBF = div(Tf,DB,0,1);
Jij = -div(b,DBF,1,1);

DB = div(Tnu0,D,0,1); DBF = div(Tf,DB,0,1);
Eij = -div(b,DBF,1,1);
end 