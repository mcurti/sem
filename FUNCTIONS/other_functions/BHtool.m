% *** Load material properties (BH-curve) from a text file, that can be
%     imported manually from FLUX/Altair Material Manager.
%     Copy/Paste/Adapt l(11-12) accordingly to add new material options.
% *** Takes in input the Magnetic field modulus B=sqrt(Bx^2+By^2), 
%     and the non-linear material BH curve choice, 1,2,3...
% *** Return the remanent field source Br and the relative permeability mu_r
%     Copyright L.A.J. Friedrich    24/10/2017
function [Br, mur] = BHtool (Bmod,option)
% Material selection
if      option==1
    [H,B] = textscan('BH_Cogent_M800_50A_50Hz.txt','%f');Bsat=(B(end-1)+B(end))./2; % Saturation field, no tangent further
%     Bmin=0.09379;               % If lower tha Bmin, interpolation gives negative Br
    Bmin=0.1379;
%     Interpolation
    mu0 = 4*pi*1e-7;
    SB  = spline(B,H);
    x0  = Bmod;
    x0(x0>Bsat) = Bsat;
    x0(x0<Bmin) = Bmin;
elseif  option==2
    fileID = fopen('BH_Cogent_M270_35A_50Hz.txt');
    out = textscan(fileID,'%f %f');  H = out{1}; B = out{2};
    fclose(fileID);
%     B(end) = (H(end)-H(end-1))*4*pi*1e-7 + B(end-1);
%     Bsat= B(end);%(B(end-1)+B(end))./2; % Saturation field, no tangent further
%     Bmin=0.09379;               % If lower tha Bmin, interpolation gives negative Br
%     Bmin=0.705;
%     Interpolation
%     mur = 1.6e4; 
    mu0 = 4*pi*1e-7;
%     Js = 2.1;
%     bh = @(H) mu0*H+2*Js/pi.*atan((pi*(mur-1)*mu0*H)./(2*Js));
    ik = B>0 & B<0.9;
    B(ik) = []; H(ik) = [];
%     B = bh(H);
    SB   = spline(B,[1/(1.6e4*mu0); H; 1/mu0]);%
    SBex = fnxtr(SB,2);
    x0  = Bmod;
%     x0(x0>Bsat) = Bsat;
%     x0(x0<Bmin) = Bmin;
else
    warning ('Undefined material option - Please import new material data as explained in description');
end


% f   = @(x) ppval(SBex,x);
fp  = fnder(SBex,1);
% fp2  = fnder(SBex,2);
% fp3  = fnder(SBex,3);

% Bf = B; Hf = H;
% FP2 = ppval(fp2,B);
% FP3 = ppval(fp3,B); FP3p = FP3(1:10); FP3p(abs(FP3p)>1e4) = 2e8; FP3(1:10) = FP3p;
% ki1 =abs(FP3)>1.5e8; %ki2 = FP2<100; ki = or(ki1 , ki2); %ki(end) = false;
% Bf(ki) = []; Hf(ki) = [];
% SBf = spline(Bf,[1/(2e4*mu0); Hf; 1/mu0]);
% SBfex = fnxtr(SBf,2);
% fpf  = fnder(SBfex,1);
f   = @(x) ppval(SBex,x);
df  = ppval(fp,x0); 
% Output
mur = 1./df./mu0;     
Br  = x0 - f(x0)./df; 
end