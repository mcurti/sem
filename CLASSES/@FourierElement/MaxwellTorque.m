% Function to evaluate the forces
function Fx = MaxwellTorque(obj)
if isempty(obj.c1)
    error('The coefficients for the solution are not uploaded!')
end

% Intialisation
Q          = obj.Edata.Harmonics;
type       = obj.Edata.type;
y1 = obj.Edata.heights(2); y2 = obj.Edata.heights(1);
C1 = obj.c1(1:Q);   C2 = obj.c2(1:Q);
C3 = obj.c3(1:Q);   C4 = obj.c4(1:Q);
C0 = obj.s(1:Q);
%==============================================================

switch type
    case 'cartesian'
        
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
        
    case 'polar'
        
        y_j = linspace(y1,y2,10);
        Fx  = zeros(1,10);
        
        for k = 1:10
        [hp, hn, bp, bn] = hbnp_r(obj.w_n,y_j(k),obj.ys(1),obj.ys(2));
                
        br =  ((bp.*C1'.*obj.w_n) + (bn.*C2'.*obj.w_n))/y_j(k); 
        ar = -((bp.*C3'.*obj.w_n) + (bn.*C4'.*obj.w_n))/y_j(k);
        
        at = -((hp.*C1') + (hn.*C2'));
        bt = -((hp.*C3') + (hn.*C4'));
        
        Fx(k) = (sum(ar.*at) +  sum(br.*bt))/(pi*4e-7)*(y_j(k))^2*(obj.Edata.tau);
        end
        Fx = mean(Fx);
end

end