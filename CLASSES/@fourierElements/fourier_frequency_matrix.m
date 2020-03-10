%==========================================================================
% fourier_frequency_matrix.m
% Created: 09.02.2018 - 13:26:17
% By: M. Curti
%
% Function to build the part of the Global matrix which coresponds to
% Fourier sapce. For this function the ProblemData structure is needed at 
% the input
%==========================================================================
function Efrequency = fourier_frequency_matrix...
                                     (obj,ProblemData,Geometry,FourierData)
% initiating the general data
conn_lines = FourierData.connectedLines;
conn_El    = FourierData.connectedElements;
metrics    = Geometry.metrics;
conn_point = Geometry.lines.connectivity(conn_lines,:);
x_fourier  = cell(1,numel(conn_lines));
w_fourier  = cell(1,numel(conn_lines));
Nel        = numel(conn_lines);% Number of elements
Nnodes     = zeros(1,Nel);     % Number of nodes in each element
Nfel       = obj.Edata.Nfel;   % Number of Fourier Elements
Q          = obj.Edata.Harmonics;
hf         = obj.Edata.heights;% Heights of Fourier Elements
FourierSize= Nfel*Q*4 + 1;
type       = obj.Edata.type;
%===============================================================
for k = 1:Nel
    x = Geometry.lines.vector{conn_lines(k)}(1,:);
    y = Geometry.lines.vector{conn_lines(k)}(2,:);
    
    theta = atan2(y,x);
    
    theta(theta< -pi/2) = 2*pi - abs(theta(theta< -pi/2));
    
    
    if strcmp(type,'polar')
        xi = theta;
    else
        xi = x;
    end
    
    x_fourier{k}  = xi;
    Nnodes(k)     = numel(x_fourier{k});
    l = conn_El(k);
    %------------------------------------------------------
    neul       = int_el(metrics.Neu2{l},1,'lines');
    neup       = int_el(metrics.Neu2{l},[1 2],'points');
    %------------------------------------------------------
    neu2       = [neup(1) neul{1} neup(2)];
    w_fourier{k}  = ...
        Geometry.xi.w_for_all_lines{conn_lines(k)}'.*neu2;
end
NnodesSum = sum(Nnodes);  % Total number of nodes;
% SEM to fourier transformation
[Ic0, Ia, Ib]   = lagrange2fourier_int(x_fourier,w_fourier, Q);

Ic0(1) = Ic0(1)+Ic0(end);
Ia(:,1) = Ia(:,1)+Ia(:,end);
Ib(:,1) = Ib(:,1)+Ib(:,end);
% Initiating the Efrequency matrix
unknowns       = ProblemData.inElementsUnknowns + ...
    ProblemData.linesUnknowns + ProblemData.pointsUnknowns; 
    % The unknowns in the SEM model
    
% Fourier size
Efrequency = zeros(FourierSize,unknowns);

sin_index       = (1:2*Q) + unknowns;
cos_index       = (1:2*Q) + 2*Q + unknowns;
dc_index        = [1 2] + unknowns + Q*4*Nfel;

unique_points = unique(conn_point,'stable');
unique_points(end) = [];
index_space = [ProblemData.lineVectorLocation{conn_lines}...
                           ProblemData.pointVectorLocation(unique_points)];

% Index for the Efrequency lines and points
P = zeros(1,Nel+1); f_lines = zeros(1,NnodesSum-Nel*2);
tmp_inc = 1:NnodesSum-Nel+1;

currentIndex = 0;
lineIndex    = 0;
for k = 1:Nel+1
    P(k) = tmp_inc(1 + currentIndex);
    if k <= Nel
        f_lines((1:Nnodes(k)-2) + lineIndex) = tmp_inc((2:Nnodes(k)-1) ...
                                                           + currentIndex);
                                                       
        currentIndex = currentIndex + Nnodes(k)-1;
        lineIndex    = lineIndex + Nnodes(k)-2;
    end
end

index_fourier = [f_lines P(1:end-1)];


% Building the boundary condition matrix


FourierMatrix = zeros(FourierSize);
K = @(k) k*4*Q;
for k = 1:Nfel
    if strcmp(type,'cartesian')
        [p_y1, n_y1] = ...
                           np_exp(obj.w_n,hf(k,1),obj.ys(k,1),obj.ys(k,2));
        [p_y2, n_y2] = ...
                           np_exp(obj.w_n,hf(k,2),obj.ys(k,1),obj.ys(k,2));
        c = hf(1,1);
    else
        [~, ~, p_y1, n_y1] = ...
                           hbnp_r(obj.w_n,hf(k,1),obj.ys(k,1),obj.ys(k,2));
        [p_y2, n_y2, ~, ~] = ...
                           hbnp_r(obj.w_n,hf(k,2),obj.ys(k,1),obj.ys(k,2));
        c = 0;
                       
    end
    
    if k == 1
        FourierMatrix(1:Q,1:2*Q) = -[diag(p_y1)  diag(n_y1)];
        FourierMatrix((1:Q) + 2*Q,(1:2*Q) + 2*Q) = ...
                                                  -[diag(p_y1) diag(n_y1)];
                                                       
        % The DC component with SEM
        %                    FourierMatrix(Nfel*Q*4 + 1,Nfel*Q*4 + 1) = ...
        %                                                         -hf(k,1);
    elseif k == 2
        %Az
        %sine
        FourierMatrix((1:Q)+Q,(1:6*Q)) = [diag(p_y1), ...
           diag(n_y1), zeros(Q,2*Q), -diag(p_y1), -diag(n_y1)];
        %cosine
        FourierMatrix((1:Q)+3*Q,(1:6*Q)+2*Q) = [diag(p_y1) ...
           diag(n_y1), zeros(Q,2*Q), -diag(p_y1), -diag(n_y1)];
        
        %Hx
        %sine
        FourierMatrix((1:Q)+4*Q,(1:6*Q)) = [diag(p_y1), ...
          -diag(n_y1), zeros(Q,2*Q), -diag(p_y1),  diag(n_y1)];
        
        %cosine
        FourierMatrix((1:Q)+6*Q,(1:6*Q)+2*Q) = [diag(p_y1), ...
           -diag(n_y1), zeros(Q,2*Q), -diag(p_y1), diag(n_y1)];
        % The DC components to space
        %          FourierMatrix(end-1,[K(k)+2 end]) =  [hf(k,1) -hf(k,1)];
        
    end
    
    
    
    if k == Nfel
        FourierMatrix(K(k) - 2*Q + (-Q+1:0), K(k) - 2*Q + (-2*Q+1:0)) = ...
                                                  [diag(p_y2) -diag(n_y2)];
        
        FourierMatrix(K(k) + (-Q+1:0),K(k) + (-2*Q+1:0)) = ...
                                                  [diag(p_y2) -diag(n_y2)];
        % The DC components to space
        FourierMatrix(end,K(k) + 1) = -1*c;
    end
    
end
% Building the Efrequency based on boundary conditions

% Sine term for the shared line
% SEM side
Efrequency(sin_index(1:end/2)-unknowns,index_space) = Ia(:,index_fourier);

% Cosine terms
% SEM side
Efrequency(cos_index(1:end/2)-unknowns,index_space) = Ib(:,index_fourier);

% DC components
Efrequency(dc_index(1)-unknowns,index_space)        = Ic0(index_fourier);
Efrequency = [Efrequency FourierMatrix];

end