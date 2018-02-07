% *** Load material properties (BH-curve) from a text file, that can be
%     imported manually from FLUX/Altair Material Manager.
%     Copy/Paste/Adapt l(11-12) accordingly to add new material options.
% *** Takes in input the Magnetic field modulus B=sqrt(Bx^2+By^2), 
%     and the non-linear material BH curve choice, 1,2,3...
% *** Return the remanent field source Br and the relative permeability mu_r
%     Copyright L.A.J. Friedrich    24/10/2017
function [Br, mur_inc] = MUBrtool (mur,option)
% Material selection
if      option==1
    [H,B] = textread('BH_Cogent_M800_50A_50Hz.txt');
elseif  option==2
    [H,B] = textread('BH_Cogent_M270_35A_50Hz.txt'); 
else 
    warning ('Undefined material option - Please import new material data as explained in description');
end
Bsat=(B(end-1)+B(end))./2; % Saturation field, no tangent further
Bmin=0.09379;               % If lower tha Bmin, interpolation gives negative Br
% Interpolation
mu0 = 4*pi*1e-7;   
SB  = spline(B,H);
x0  = Bmod;
x0(x0>Bsat) = Bsat;
x0(x0<Bmin) = Bmin;

f   = @(x) ppval(SB,x);
df  = derivest(f,x0); 
% Output
mur_inc = 1./df./mu0;     
Br  = -f(x0)./df + x0; 
end