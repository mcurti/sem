function out = int_el(varargin)
%------------------------------------------------------------------
% Interior nodes values of an element with data
% intElement - returns the interior elements in El
% lines      - returns the lines in El
% points     - returns the points in El
% adress     - returns the addresses of interior elements, lines and points
% from El;
%------------------------------------------------------------------
    mode = varargin{nargin};
    El   = varargin{1};
    [M, N] = size(El);
    switch mode
        case 'intElement'
            out = El(2:M-1,2:N-1);
        case 'lines'
            lines = varargin{2};
            out = cell(1,length(lines));
            for k = 1:length(lines)
                if lines(k) == 1
                    out{k} = reshape(El(1,2:end-1),1,N-2);
                elseif lines(k) == 2
                    out{k} = reshape(El(2:end-1,end),1,M-2);
                elseif lines(k) == 3
                    out{k} = reshape(El(end,2:end-1),1,N-2);
                elseif lines(k) == 4
                    out{k} = reshape(El(2:end-1,1),1,M-2);
                end
            end
        case 'points'
            points = varargin{2};
            out = zeros(1,length(points));
            for k = 1:length(points)
                if points(k) == 1
                    out(k) = El(1,1);
                elseif points(k) == 2
                    out(k) = El(1,end);
                elseif points(k) == 3
                    out(k) = El(end,end);
                elseif points(k) == 4
                    out(k) = El(end,1);
                end
            end
        case 'address'
            [xx, yy] = meshgrid(1:N, 1:M);
             
             yel  = find(abs(xx)<N & abs(yy)<M & abs(xx)>1 & abs(yy)>1);
             
             yc     = [find(xx==1 & yy==1) find(xx==N & yy==1)...
                       find(xx==N & yy==M) find(xx==1 & yy==M)];
             yb{1}   = find(xx>1 & xx<N & yy==1)';
             yb{2}   = find(xx==N & yy<M & yy>1)';
             yb{3}   = find(xx>1 & xx<N & yy==M)';  
             yb{4}   = find(xx==1 & yy<M & yy>1)';
             out{1} = yel; out{2} = yb; out{3} = yc;
        
        case 'lines_address'
            [xx, yy] = meshgrid(1:N, 1:M);
            
             yb{1}   = find(xx>=1 & xx<=N & yy==1)';
             yb{2}   = find(xx==N & yy<=M & yy>=1)';
             yb{3}   = find(xx>=1 & xx<=N & yy==M)';  
             yb{4}   = find(xx==1 & yy<=M & yy>=1)';
             out = yb;
    end 
end
