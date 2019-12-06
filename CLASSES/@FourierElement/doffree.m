% Function to evaluate the forces
function out = doffree(obj)


% Intialisation
Q          = obj.Edata.Harmonics;
type       = obj.Edata.type;
y1 = obj.Edata.heights(2); y2 = obj.Edata.heights(1);

mover_offset  = obj.Edata.x_start*2*pi/obj.Edata.tau;
alpha         = mover_offset*pi/180*(1:Q);

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
        
        out = (sum(ax.*by) -  sum(bx.*ay))/(pi*4e-7)*100/obj.w_n(1)*pi;
        
    case 'polar'
        
        % y1 - stator
        % y2 - rotor
        
       
        Ty   = diag((obj.ys(2)/y2).^obj.w_n./...
            ((y1/y2*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (y2/y1*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tm   = diag((obj.ys(2)/y1).^obj.w_n./...
            ((y1/y2*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (y2/y1*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tg   = diag((y1/obj.ys(1)).^obj.w_n./...
            ((y1/y2*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (y2/y1*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tp   = diag((y2/obj.ys(1)).^obj.w_n./...
            ((y1/y2*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (y2/y1*obj.ys(2)/obj.ys(1)).^obj.w_n));
        
        C1t = Ty*azss - Tm*azsr;
        C2t = Tg*azsr - Tp*azss;
        C3t = Ty*azcs - Tm*azcr;
        C4t = Tg*azcr - Tp*azcs;
        
        out = [azsr; azss; azcr; azcs];
        
end

end