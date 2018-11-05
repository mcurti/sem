%==========================================================================
% SEM.m
% Created: 12.02.2018 - 08:23:21
% By: M. Curti
%
% Function to evaluate the spatial values of the solution within
% the defined domain
%==========================================================================
function [x_grid, y_grid, f_solution] = fourier2space(obj, dx, dy, El)

if isempty(obj.c1)
    error('The coefficients for the solution are not upladed!')
end

% Intialisation
type = obj.Edata.type;
Q  = obj.Edata.Harmonics; 
x1 = 0;       x2 = x1 + obj.Edata.tau;
y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);

C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
C0 = obj.s((1:2*Q) + (El-1)*2*Q);
%==============================================================

x_i = linspace(x1, x2, dx); y_j = linspace(y1, y2, dy);
% theta = linspace(0, abs(obj.Edata.tau), dx);

[x_grid, y_grid] = meshgrid(x_i,y_j);
x_grid = fliplr(x_grid);
if strcmp(type,'cartesian')
    a  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1) + ...
        exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2) - ...
        repmat(C0((1:Q)+Q)./obj.w_n',dy,1);
    b  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3) + ...
        exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4) + ...
        repmat(C0(1:Q)./obj.w_n',dy,1);
    
    c0 = obj.Az0 + obj.Bx0*y_grid;
else
    rp = zeros(dy,Q); rn = rp;
    for k = 1:dy
        for h = 1:Q
            rp(k,h) = (y_j(k)/obj.ys(El,1))^obj.w_n(h);
            rn(k,h) = (y_j(k)/obj.ys(El,2))^(-obj.w_n(h));
        end
    end
    a = (rp*diag(C1) + rn*diag(C2)); b = (rp*diag(C3) + rn*diag(C4));
    c0 = obj.Bx0*log(y_grid) + obj.Az0;
end


f_solution = c0 + a*sin(obj.w_n*x_i) + b*cos(obj.w_n*x_i);
end
