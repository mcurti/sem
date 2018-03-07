%==========================================================================
% fourier_frequency_matrix.m
% Created: 09.02.2018 - 13:26:17
% By: M. Curti
%
% Function to build the part of the Global matrix which coresponds to
% Fourier sapce. For this function the ProblemData structure is needed at
% the input
%==========================================================================
function Efrequency = fourier_frequency_matrix(obj,ProblemData)
% initiating the general data
conn_lines{1} = obj.SEMdata.Lines_top;
conn_lines{2} = obj.SEMdata.Lines_bot;
% metrics       = obj.SEMdata.metrics;
% conn_El{1}    = obj.SEMdata.Elements_top;
% conn_El{2}    = obj.SEMdata.Elements_bot;
conn_point{1} = obj.SEMdata.lines.connectivity(conn_lines{1},:);
conn_point{2} = obj.SEMdata.lines.connectivity(conn_lines{2},:);
Nel(1)        = numel(conn_lines{1}); % Number of elements
Nel(2)        = numel(conn_lines{2}); % Number of elements
Nnodes{1}     = zeros(1,Nel(1));      % Number of nodes in each element
Nnodes{2}     = zeros(1,Nel(2));      % Number of nodes in each element
NnodesSum     = zeros(1,2);
% Nfel       = obj.Edata.Nfel;   % Number of Fourier Elements
Q          = obj.Edata.Harmonics;
hf         = obj.Edata.heights;% Heights of Fourier Elements
FourierSize= Q*4 + 2;
type       = obj.Edata.type;
unknowns       = ProblemData.inElementsUnknowns + ...
    ProblemData.linesUnknowns + ProblemData.pointsUnknowns;
switch type
    case 'cartesian'
        [p_y1, n_y1] = ...
            np_exp(obj.w_n,hf(1),obj.ys(1),obj.ys(2));
        [p_y2, n_y2] = ...
            np_exp(obj.w_n,hf(2),obj.ys(1),obj.ys(2));
        c = [1 2];
%         c1 = hf(2,1);
    case 'polar'
        [~, ~, p_y1, n_y1] = ...
            hbnp_r(obj.w_n,hf(1),obj.ys(1),obj.ys(2));
        [~, ~, p_y2, n_y2] = ...
            hbnp_r(obj.w_n,hf(2),obj.ys(1),obj.ys(2));
        c = [1 2];
%         c1 = [log(hf(1)) 1];
end
sin_index       = (1:2*Q) + unknowns;
cos_index       = (1:2*Q) + 2*Q + unknowns;
dc_index        = c + unknowns + Q*4;

% Initiating the Efrequency matrix
% The unknowns in the SEM model
% Building the boundary condition matrix


FourierMatrix = zeros(FourierSize);

FourierMatrix(1:Q,1:2*Q)                       =  -[diag(p_y1) diag(n_y1)]*hf(1);
FourierMatrix((1:Q) + 2*Q,(1:2*Q) + 2*Q)       =  -[diag(p_y1) diag(n_y1)]*hf(1);


FourierMatrix(2*Q + (-Q+1:0),2*Q + (-2*Q+1:0)) =  -[diag(p_y2) diag(n_y2)]*hf(2);
FourierMatrix(4*Q + (-Q+1:0),4*Q + (-2*Q+1:0)) =  -[diag(p_y2) diag(n_y2)]*hf(2);

% The DC components to space
FourierMatrix(FourierSize-1,FourierSize-[1 0]) =  -[log(hf(1)) 1];
FourierMatrix(FourierSize,FourierSize-[1 0])   =  -[log(hf(2)) 1];
% Fourier size
Efrequency = zeros(FourierSize,unknowns);
%===============================================================

for q = 1:2
    xi_fourier  = cell(1,Nel(q));
    w_fourier   = cell(1,Nel(q));
%     nl = 1+(2*q-2);
    for k = 1:Nel(q)
        x = obj.SEMdata.lines.vector{conn_lines{q}(k)}(1,:);
        y = obj.SEMdata.lines.vector{conn_lines{q}(k)}(2,:);
        
        theta = atan2(y,x);
        
        theta(theta< -pi/2) = 2*pi - abs(theta(theta< -pi/2));
        
        
        if strcmp(type,'polar')
            if k == 1
            theta_start = theta(1);
            theta_sign  = sign(theta(end) - theta(1));
            end
            xi = (theta - theta_start)*theta_sign;
        else
            xi = x;
        end
        
        xi_fourier{k}  = xi;
        Nnodes{q}(k)     = numel(xi_fourier{k});
        
        w_fourier{k}  = ...
            obj.SEMdata.xi.w_for_all_lines{conn_lines{q}(k)}'*((max(xi) - min(xi))/2);  %
    end
    NnodesSum(q) = sum(Nnodes{q});  % Total number of nodes;
    % SEM to fourier transformation
    [Ic0, Ia, Ib]   = lagrange2fourier_int(xi_fourier,w_fourier, Q);
    
    Ic0([1 end])  = Ic0(1)+Ic0(end);
    Ia(:,[1 end]) = repmat(Ia(:,1)+Ia(:,end),1,2);
    Ib(:,[1 end]) = repmat(Ib(:,1)+Ib(:,end),1,2);
    
    unique_points = unique(conn_point{q},'stable');
    unique_points(end) = [];
    index_space = [ProblemData.lineVectorLocation{conn_lines{q}}...
        ProblemData.pointVectorLocation(unique_points)];
    
    % Index for the Efrequency lines and points
    P = zeros(1,Nel(q)+1); f_lines = zeros(1,NnodesSum(q)-Nel(q)*2);
    tmp_inc = 1:NnodesSum(q)-Nel(q)+1;
    
    currentIndex = 0;
    lineIndex    = 0;
    for k = 1:Nel(q)+1
        P(k) = tmp_inc(1 + currentIndex);
        if k <= Nel(q)
            f_lines((1:Nnodes{q}(k)-2) + lineIndex) =...
                tmp_inc((2:Nnodes{q}(k)-1) + currentIndex);
            currentIndex = currentIndex + Nnodes{q}(k)-1;
            lineIndex    = lineIndex + Nnodes{q}(k)-2;
        end
    end
    
    index_fourier = [f_lines P(1:end-1)];
    
    if q == 1; dc = 1; else; dc = 1; end
    
    % Building the Efrequency based on boundary conditions
    
    % Sine term for the shared line
    % SEM side
    Efrequency(sin_index((1:Q)+Q*(q-1))-unknowns,index_space) =  Ia(:,index_fourier)*dc;
    
    % Cosine terms
    % SEM side
    Efrequency(cos_index((1:Q)+Q*(q-1))-unknowns,index_space) =  Ib(:,index_fourier)*dc;
    
    % DC components
    Efrequency(dc_index(q)-unknowns,index_space)              = Ic0(index_fourier)*dc;
end
Efrequency = [Efrequency, FourierMatrix];
end