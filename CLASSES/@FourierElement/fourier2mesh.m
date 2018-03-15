function f_solution = fourier2mesh(obj, x, y,El)

if isempty(obj.c1)
    error('The coefficients for the solution are not upladed!')
end
% Intialisation

type = obj.Edata.type;
Q    = obj.Edata.Harmonics;
C1   = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
C3   = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
C0   = obj.s((1:2*Q) + (El-1)*2*Q);
%==============================================================
points_number = numel(x);

f_solution = zeros(points_number,1);

if strcmp(type,'cartesian')
    for k = 1:points_number
        a  = exp( obj.w_n*(y(k)-obj.ys(El,1))).*C1' + ...
            exp(-obj.w_n*(y(k)-obj.ys(El,2))).*C2'+ C0((1:Q)+Q)'./obj.w_n;
        b  = exp( obj.w_n*(y(k)-obj.ys(El,1))).*C3' + ...
            exp(-obj.w_n*(y(k)-obj.ys(El,2))).*C4'+ C0(1:Q)'./obj.w_n;
        c0 = obj.Az0 + obj.Bx0*y(k);
        
        f_solution(k) = c0 + sum(a.*sin(obj.w_n*x(k)) + ...
            b.*cos(obj.w_n*x(k)));
    end
else
%     y = pi-y*pi/180;
    
    [t, r] = cart2pol(x,y);
    t = pi-t;
    for k = 1:points_number
        
        [~, ~, p_y1, n_y1] = ...
                            hbnp_r(obj.w_n,r(k),obj.ys(El,1),obj.ys(El,2));
        a  = p_y1.*C1' + n_y1.*C2';
        b  = p_y1.*C3' + n_y1.*C4';
        c0 = obj.Az0 + obj.Bx0*log(r(k));
        
        f_solution(k) = c0 + sum(a.*sin(obj.w_n*t(k)) + ...
            b.*cos(obj.w_n*t(k)));
    end
end
end
