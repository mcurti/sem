%==========================================================================
% fourier_space_matrix.m
% Created: 09.02.2018 - 10:57:54
% By: M. Curti
%
% Function to build the part of the Global matrix which coresponds to
% Fourier evaluated in spatial domain. For this function the ProblemData
% structure is needed at the input
%==========================================================================
function [Espace, interface_line] = fourier_space_matrix(obj,...
    ProblemData,Geometry,FourierData)
% Initiating the general data
conn_lines    = FourierData.connectedLines;
metrics       = Geometry.metrics;
conn_El       = FourierData.connectedElements;
conn_point    = Geometry.lines.connectivity(conn_lines,:);
Nel           = numel(conn_lines); % Number of elements
Nnodes        = zeros(1,Nel);      % Number of nodes in each element
Nfel          = obj.Edata.Nfel;    % Number of Fourier Elements
Q             = obj.Edata.Harmonics;
type          = obj.Edata.type;

for k = 1:Nel
    Nnodes(k) = numel(Geometry.lines.vector{conn_lines(k)}(1,:));
end
NnodesSum = sum(Nnodes);  % Total number of nodes;

% Building the quadrature and abscisas vector

xi_fourier     = zeros(1,NnodesSum-Nel+1);
w_fourier     = zeros(1,NnodesSum-Nel+1);

currentIndex = Nnodes(1);
for k = 1:Nel
    %------------------------------------------------------
    neul       = int_el(metrics.Neu2{conn_El(k)},1,'lines');
    neup       = int_el(metrics.Neu2{conn_El(k)},[1 2],'points');
    neu2       = [neup(1) neul{1} neup(2)];
    l = conn_lines(k);
    x = Geometry.lines.vector{l}(1,:);
    y = Geometry.lines.vector{l}(2,:);
    
    theta = atan2(y,x);
    
    theta(theta< -pi/2) = 2*pi - abs(theta(theta< -pi/2));
    
    
    if strcmp(type,'polar')
        xi = theta;
    else
        xi = x;
    end
    
    if k==1
        xi_fourier(1:Nnodes(1)) = xi;
        w_fourier(1:Nnodes(1)) = Geometry.xi.w_for_all_lines{l}'.*neu2;
    else
        xi_fourier((0:Nnodes(k)-1) + currentIndex) = xi;
        w_fourier((0:Nnodes(k)-1) + currentIndex) = ...
            w_fourier((0:Nnodes(k)-1) + currentIndex) + ...
            Geometry.xi.w_for_all_lines{l}'.*neu2;
        
        currentIndex = currentIndex + Nnodes(k) - 1;
    end
end
w_fourier([1 end]) = w_fourier(1) + w_fourier(end);
% Removing the last entry in the vectors
%             x_fourier(end) = []; w_fourier(end) = [];

% Initiating the Espace matrix
unknowns = ProblemData.inElementsUnknowns + ProblemData.linesUnknowns + ...
    ProblemData.pointsUnknowns;
Espace   = zeros(unknowns,Q*4*Nfel + 1);

% Initiating the entries, the exponentials are necessary only
% for the first layer, since other layers are not connected to
% the SEM part


if strcmp(type,'cartesian')
    [p_exp, n_exp] = np_exp(obj.w_n,obj.Edata.heights(1,1),...
        obj.ys(1,1),obj.ys(1,2));
    % Computing the block matrix in the Espace matrix
    Espace_block = - diag(w_fourier(1:end-1))*...
        [sin(xi_fourier(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
       - sin(xi_fourier(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
         cos(xi_fourier(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
       - cos(xi_fourier(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
        zeros(numel(xi_fourier(1:end-1)),Q*4*(Nfel-1)),...
        ones(numel(xi_fourier(1:end-1)),1)];
else
    [h_p, h_n, ~, ~] = hbnp_r(obj.w_n,obj.Edata.heights(1,1),...
        obj.ys(1,1),obj.ys(1,2));
    % Computing the block matrix in the Espace matrix
    Espace_block = -(diag(w_fourier(1:end-1)))/obj.Edata.heights(1,1)*...
        [  sin(xi_fourier(1:end-1)'*obj.w_n')*diag(h_p)*diag(obj.w_n), ...
         - sin(xi_fourier(1:end-1)'*obj.w_n')*diag(h_n)*diag(obj.w_n), ...
           cos(xi_fourier(1:end-1)'*obj.w_n')*diag(h_p)*diag(obj.w_n), ...
         - cos(xi_fourier(1:end-1)'*obj.w_n')*diag(h_n)*diag(obj.w_n), ...
        zeros(numel(xi_fourier(1:end-1)),Q*4*(Nfel-1)),...
        zeros(numel(xi_fourier(1:end-1)),1)];
end






% Computing the Espace_block address in the Espace
unique_points = unique(conn_point,'stable');
unique_points(end) = [];

index_space = [ProblemData.lineVectorLocation{conn_lines}...
    ProblemData.pointVectorLocation(unique_points)];

% Index for the Espace lines and points
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
% Moving the block to the matrix
Espace(index_space,:) = Espace_block(index_fourier,:);


interface_line.x = xi_fourier;
interface_line.w = w_fourier;
end

