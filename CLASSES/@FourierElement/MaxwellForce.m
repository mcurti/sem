% Function to evaluate the forces
function Fx = MaxwellForce(obj)
if isempty(obj.c1)
    error('The coefficients for the solution are not uploaded!')
end

% Intialisation
Q  = obj.Edata.Harmonics;

y1 = obj.Edata.heights(1,2); y2 = obj.Edata.heights(1,1);
C1 = obj.c1(1:Q);   C2 = obj.c2(1:Q);
C3 = obj.c3(1:Q);   C4 = obj.c4(1:Q);
C0 = obj.s(1:Q);
%==============================================================           x_i = linspace(x1, x2, dx);
y_j = y1 + (y2-y1)/2;


ax  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') - ...
    exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n');
bx  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') - ...
    exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n');


ay  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') + ...
    exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n');
by  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') + ...
    exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n') - ...
    C0.*obj.w_n';

Fx = (sum(ax.*by) -  sum(bx.*ay))/(pi*4e-7)*100/obj.w_n(1)*pi;

end