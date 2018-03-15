%==========================================================================
% fourier_space_matrix.m
% Created: 14.02.2018 - 16:27:10
% By: M. Curti
%
% Function to build the part of the Global matrix which coresponds to
% Fourier evaluated in spatial domain. For this function the ProblemData
% structure is needed at the input
%==========================================================================
function Espace = fourier_space_matrix(obj,ProblemData)
% Initiating the general data
conn_lines{1} = obj.SEMdata.Lines_top;
conn_lines{2} = obj.SEMdata.Lines_bot;
metrics       = obj.SEMdata.metrics;
conn_El{1}    = obj.SEMdata.Elements_top;
conn_El{2}    = obj.SEMdata.Elements_bot;
type          = obj.Edata.type;
conn_point{1} = obj.SEMdata.lines.connectivity(conn_lines{1},:);
conn_point{2} = obj.SEMdata.lines.connectivity(conn_lines{2},:);
Nel(1)        = numel(conn_lines{1}); % Number of elements
Nel(2)        = numel(conn_lines{2}); % Number of elements
Nnodes{1}     = zeros(1,Nel(1));      % Number of nodes in each element
Nnodes{2}     = zeros(1,Nel(2));      % Number of nodes in each element
hf            = obj.Edata.heights;
Q             = obj.Edata.Harmonics;

for k = 1:Nel(1)
    Nnodes{1}(k) = numel(obj.SEMdata.lines.vector{conn_lines{1}(k)}(1,:));
end
for k = 1:Nel(2)
    Nnodes{2}(k) = numel(obj.SEMdata.lines.vector{conn_lines{2}(k)}(1,:));
end
NnodesSum{1} = sum(Nnodes{1});  % Total number of nodes;
NnodesSum{2} = sum(Nnodes{2});  % Total number of nodes;

% Building the quadrature and abscisas vector

xi_fourier{1}     = zeros(1,NnodesSum{1}-Nel(1)+1);
xi_fourier{2}     = zeros(1,NnodesSum{2}-Nel(2)+1);
w_fourier{1}      = zeros(1,NnodesSum{1}-Nel(1)+1);
w_fourier{2}      = zeros(1,NnodesSum{2}-Nel(2)+1);


for q = 1:2
    currentIndex = Nnodes{q}(1);
    for k = 1:Nel(q)
        %------------------------------------------------------
        
        nl = 1+(2*q-2);
%         nl = 3;
        if q == 1
            sn = -1;
        else
            sn = 1;
        end
        lines_index = int_el(metrics.Neu2{conn_El{q}(k)},nl,'lines_address');
        neu2      = metrics.Neu2{conn_El{q}(k)}(lines_index{nl})*sn;
        l = conn_lines{q}(k);
        x = obj.SEMdata.lines.vector{l}(1,:);
        y = obj.SEMdata.lines.vector{l}(2,:);
        
        theta = atan2(y,x);
        
        theta(theta< -pi/2) = 2*pi - abs(theta(theta< -pi/2));
        
        
        if strcmp(type,'polar')
            xi = theta;
        else
            xi = x;
        end
        
        if k==1
            xi_fourier{q}(1:Nnodes{q}(1)) = xi;
            w_fourier{q}(1:Nnodes{q}(1)) = obj.SEMdata.xi.w_for_all_lines{l}'.*neu2;
        else
            xi_fourier{q}((0:Nnodes{q}(k)-1) + currentIndex) = xi;
            w_fourier{q}((0:Nnodes{q}(k)-1) + currentIndex) = ...
                w_fourier{q}((0:Nnodes{q}(k)-1) + currentIndex) + ...
                obj.SEMdata.xi.w_for_all_lines{l}'.*neu2;
            
            currentIndex = currentIndex + Nnodes{q}(k) - 1;
        end
    end
    
    if strcmp(type,'polar')
        theta_start = xi_fourier{q}(1);
        theta_sign  = sign(xi_fourier{q}(end) - xi_fourier{q}(1));
        xi_fourier{q} = (xi_fourier{q} - theta_start)*theta_sign;
    end
    w_fourier{q}([1 end]) = w_fourier{q}(1) + w_fourier{q}(end);
end
% Removing the last entry in the vectors
%             x_fourier(end) = []; w_fourier(end) = [];

% Initiating the Espace matrix

unknowns = ProblemData.inElementsUnknowns + ProblemData.linesUnknowns + ...
    ProblemData.pointsUnknowns;
Espace   = zeros(unknowns,Q*4 + 2);

% Initiating the entries,
Espace_block = cell(2,1);
for k = 1:2
    if strcmp(type,'cartesian')
        [p_exp, n_exp] = np_exp(obj.w_n,obj.Edata.heights(1,1),...
            obj.ys(1,1),obj.ys(1,2));
        % Computing the block matrix in the Espace matrix
        Espace_block{k} = -diag(w_fourier{2}(1:end-1))*...
            [sin(xi_fourier{1}(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
            - sin(xi_fourier{1}(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
            cos(xi_fourier{1}(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
            - cos(xi_fourier{1}(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
            ones(numel(xi_fourier{1}(1:end-1)),1), ];
    else
        [h_p, h_n, ~, ~] = ...
            hbnp_r(obj.w_n,obj.Edata.heights(k),obj.ys(1),obj.ys(2));
    if k == 1; dc = 1; else; dc = 1; end% + 
    if k == 1; a = 0; else; a = 10*pi/180; end
        % Computing the block matrix in the Espace matrix
        Espace_block{k} = (diag(w_fourier{k}(1:end-1)))*...
            [sin(xi_fourier{k}(1:end-1)'*obj.w_n'+ a)*diag(h_p), ...
             sin(xi_fourier{k}(1:end-1)'*obj.w_n'+ a)*diag(h_n), ...
             cos(xi_fourier{k}(1:end-1)'*obj.w_n'+ a)*diag(h_p), ...
             cos(xi_fourier{k}(1:end-1)'*obj.w_n'+ a)*diag(h_n), ...
             ones(numel(xi_fourier{k}(1:end-1)),1)./hf(k)...
             , zeros(numel(xi_fourier{k}(1:end-1)),1)]*dc;
    end
    
    % Computing the Espace_block address in the Espace
    unique_points = unique(conn_point{k},'stable');
    unique_points(end) = [];
    
    index_space = [ProblemData.lineVectorLocation{conn_lines{k}}...
        ProblemData.pointVectorLocation(unique_points)];
    
    % Index for the Espace lines and points
    P = zeros(1,Nel(k)+1); f_lines = zeros(1,NnodesSum{k}-Nel(k)*2);
    tmp_inc = 1:NnodesSum{k}-Nel(k)+1;
    
    currentIndex = 0;
    lineIndex    = 0;
    for el = 1:Nel(k)+1
        P(el) = tmp_inc(1 + currentIndex);
        if el <= Nel(k)
            f_lines((1:Nnodes{k}(el)-2) + lineIndex) = tmp_inc((2:Nnodes{k}(el)-1) ...
                + currentIndex);
            currentIndex = currentIndex + Nnodes{k}(el)-1;
            lineIndex    = lineIndex + Nnodes{k}(el)-2;
        end
    end
    index_fourier = [f_lines P(1:end-1)];
    Espace(index_space,:) = Espace_block{k}(index_fourier,:);
end
end
