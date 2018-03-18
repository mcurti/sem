% Function to evaluate the forces
function Fx = MaxwellTorque(obj)
if isempty(obj.c1)
    error('The coefficients for the solution are not uploaded!')
end

% Intialisation
Q  = obj.Edata.Harmonics;

y1 = obj.Edata.heights(2); y2 = obj.Edata.heights(1);
C1 = obj.c1(1:Q);   C2 = obj.c2(1:Q);
C3 = obj.c3(1:Q);   C4 = obj.c4(1:Q);
C0 = obj.s(1:Q);
%==============================================================
y_j = y1 + (y2-y1)/2;


ax  = exp( obj.w_n*(y_j-obj.ys(1)))'*diag(C1.*obj.w_n') - ...
    exp(-obj.w_n*(y_j-obj.ys(2)))'*diag(C2.*obj.w_n');
bx  = exp( obj.w_n*(y_j-obj.ys(1)))'*diag(C3.*obj.w_n') - ...
    exp(-obj.w_n*(y_j-obj.ys(2)))'*diag(C4.*obj.w_n');


ay  = exp( obj.w_n*(y_j-obj.ys(1)))'*diag(C1.*obj.w_n') + ...
    exp(-obj.w_n*(y_j-obj.ys(2)))'*diag(C2.*obj.w_n');
by  = exp( obj.w_n*(y_j-obj.ys(1)))'*diag(C3.*obj.w_n') + ...
    exp(-obj.w_n*(y_j-obj.ys(2)))'*diag(C4.*obj.w_n') - ...
    C0.*obj.w_n';

Fx = (sum(ax.*by) -  sum(bx.*ay))/(pi*4e-7)*100/obj.w_n(1)*pi;


[h_p, h_n, rp, rn] = ...
            hbnp_r(obj.w_n,y_j,obj.ys(1),obj.ys(2));
% rp = (y_j/obj.ys(1))^obj.w_n(h);
% rn = (y_j/obj.ys(2))^(-obj.w_n(h));

ar = (rp*diag(C1) + rn*diag(C2)); br = (rp*diag(C3) + rn*diag(C4));
c0 = obj.Bx0*log(y_grid) + obj.Az0;
end