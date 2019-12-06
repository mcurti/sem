%==========================================================================
% fourier_frequency_matrix.m
% Created: 09.02.2018 - 13:26:17
% By: M. Curti
%
% Function to build the part of the Global matrix which coresponds to
% Fourier sapce. For this function the ProblemData structure is needed at
% the input
%==========================================================================
function [FourierMatrix, Efrequency] = ...
                          fourier_frequency_matrix_doffree(obj,ProblemData)
% initiating the general data
conn_lines{1} = obj.SEMdata.Lines_top;
conn_lines{2} = obj.SEMdata.Lines_bot;
metrics       = obj.SEMdata.metrics;
conn_El{1}    = obj.SEMdata.Elements_top;
conn_El{2}    = obj.SEMdata.Elements_bot;
conn_point{1} = obj.SEMdata.lines.connectivity(conn_lines{1},:);
conn_point{2} = obj.SEMdata.lines.connectivity(conn_lines{2},:);
Nel(1)        = numel(conn_lines{1}); % Number of elements
Nel(2)        = numel(conn_lines{2}); % Number of elements
Nnodes{1}     = zeros(1,Nel(1));      % Number of nodes in each element
Nnodes{2}     = zeros(1,Nel(2));      % Number of nodes in each element
NnodesSum     = zeros(1,2);
index_fourier = cell(1,2);
% Nfel       = obj.Edata.Nfel;   % Number of Fourier Elements
Q          = obj.Edata.Harmonics;
mover_offset  = obj.Edata.x_start*2*pi/obj.Edata.tau;
hf         = obj.Edata.heights;% Heights of Fourier Elements
FourierSize= Q*4 + 2;
type       = obj.Edata.type;
unknowns       = ProblemData.inElementsUnknowns + ...
    ProblemData.linesUnknowns + ProblemData.pointsUnknowns;

C12       = (1:2*Q) + unknowns;
C34       = (1:2*Q) + 2*Q + unknowns;
dc_index        = [1 2] + unknowns + Q*4;


FourierMatrix = zeros(FourierSize);

Efrequency = zeros(FourierSize,unknowns);
Espace   = zeros(unknowns,Q*4 + 2);
%===============================================================
        Ty   = diag((obj.ys(2)/hf(1)).^obj.w_n./...
            ((hf(2)/hf(1)*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (hf(1)/hf(2)*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tm   = diag((obj.ys(2)/hf(2)).^obj.w_n./...
            ((hf(2)/hf(1)*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (hf(1)/hf(2)*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tg   = diag((hf(2)/obj.ys(1)).^obj.w_n./...
            ((hf(2)/hf(1)*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (hf(1)/hf(2)*obj.ys(2)/obj.ys(1)).^obj.w_n));
        Tp   = diag((hf(1)/obj.ys(1)).^obj.w_n./...
            ((hf(2)/hf(1)*obj.ys(2)/obj.ys(1)).^obj.w_n - ...
             (hf(1)/hf(2)*obj.ys(2)/obj.ys(1)).^obj.w_n));
        ty   = 1/(log(hf(2))-log(hf(1)));
        tm   = log(hf(2))/(log(hf(2))-log(hf(1)));
%===============================================================
for k = 1:Nel(1)
    Nnodes{1}(k) = numel(obj.SEMdata.lines.vector{conn_lines{1}(k)}(1,:));
end
for k = 1:Nel(2)
    Nnodes{2}(k) = numel(obj.SEMdata.lines.vector{conn_lines{2}(k)}(1,:));
end
NnodesSum(1) = sum(Nnodes{1});  % Total number of nodes;
NnodesSum(2) = sum(Nnodes{2});  % Total number of nodes;
for q = 1:2
    xi_fourier  = cell(1,Nel(q));
    w_fourier   = cell(1,Nel(q));
    
    xi_fourier1 = zeros(1,NnodesSum(q)-Nel(q)+1);
    w_fourier1  = zeros(1,NnodesSum(q)-Nel(q)+1);
    currentIndex = Nnodes{q}(1);
    
%     nl = 1+(2*q-2);
    for k = 1:Nel(q)
        
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
        if k==1
            xi_fourier1(1:Nnodes{q}(1)) = xi;
            w_fourier1(1:Nnodes{q}(1)) = obj.SEMdata.xi.w_for_all_lines{l}'.*neu2;
        else
            xi_fourier1((0:Nnodes{q}(k)-1) + currentIndex) = xi;
            w_fourier1((0:Nnodes{q}(k)-1) + currentIndex) = ...
                w_fourier1((0:Nnodes{q}(k)-1) + currentIndex) + ...
                obj.SEMdata.xi.w_for_all_lines{l}'.*neu2;
            
            currentIndex = currentIndex + Nnodes{q}(k) - 1;
        end
        
        
    end
    
    % SEM to fourier transformation
    [Ic0, Ia, Ib]   = lagrange2fourier_int(xi_fourier, Q);
    
    Ic0([1 end])  = Ic0(1)+Ic0(end);
    Ia(:,[1 end]) = repmat(Ia(:,1)+Ia(:,end),1,2);
    Ib(:,[1 end]) = repmat(Ib(:,1)+Ib(:,end),1,2);
    
    IA{q} = Ia; IB{q} = Ib; IC0{q} = Ic0;
    unique_points = unique(conn_point{q},'stable');
    unique_points(end) = [];
    index_space{q} = [ProblemData.lineVectorLocation{conn_lines{q}}...
        ProblemData.pointVectorLocation(unique_points)];
    
%     index_space = [ProblemData.lineVectorLocation{conn_lines{q}}...
%         ProblemData.pointVectorLocation(unique_points)];
    
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
    
    index_fourier{q} = [f_lines P(1:end-1)];
    
    if q == 1
       T1 = -Tm; T2 = Tg; t1 = -ty; t2 = tm; alpha = 0;
    else
       T1 = Ty; T2 = -Tp; t1 = ty; t2 = 1-tm; alpha = -mover_offset*pi/180*(1:Q);
    end
    
    % calculate C1-4
    TC1{q}  = T1*(diag(cos(alpha))*Ia - diag(sin(alpha))*Ib);  
    TC2{q}  = T2*(diag(cos(alpha))*Ia - diag(sin(alpha))*Ib);
    TC3{q}  = T1*(diag(sin(alpha))*Ia + diag(cos(alpha))*Ib);
    TC4{q}  = T2*(diag(sin(alpha))*Ia + diag(cos(alpha))*Ib);
    TDCA{q} = t1*Ic0;
    TDCB{q} = t2*Ic0;
    

    
    %================================================================
    [h_p{q}, h_n{q}, ~, ~] = ...
            hbnp_r(obj.w_n,obj.Edata.heights(q),obj.ys(1),obj.ys(2));
    if q == 1; dc = 1; else; dc = 1; end% + 
    if q == 1; a = 0; else 
        a = (ones(numel(xi_fourier1(1:end-1)),1)*(1:Q))...
            *mover_offset*pi/180; 
    end
        % Computing the block matrix in the Espace matrix
        S{q} = diag(w_fourier1(1:end-1))*sin(xi_fourier1(1:end-1)'*obj.w_n'+ a);
        C{q} = diag(w_fourier1(1:end-1))*cos(xi_fourier1(1:end-1)'*obj.w_n'+ a);
        DC{q} = diag(w_fourier1(1:end-1))*ones(numel(xi_fourier1(1:end-1)),1)./hf(q);
end
% Building the Efrequency based on boundary conditions

% Sine term for the shared line
% SEM side
Efrequency(C12((1:Q)  )-unknowns,index_space{1}) =  TC1{1}(:,index_fourier{1});
Efrequency(C12((1:Q)+Q)-unknowns,index_space{1}) =  TC2{1}(:,index_fourier{1});
Efrequency(C12((1:Q)  )-unknowns,index_space{2}) =  TC1{2}(:,index_fourier{2});
Efrequency(C12((1:Q)+Q)-unknowns,index_space{2}) =  TC2{2}(:,index_fourier{2});

% Cosine terms
% SEM side
Efrequency(C34((1:Q)  )-unknowns,index_space{1}) =  TC3{1}(:,index_fourier{1});
Efrequency(C34((1:Q)+Q)-unknowns,index_space{1}) =  TC4{1}(:,index_fourier{1});
Efrequency(C34((1:Q)  )-unknowns,index_space{2}) =  TC3{2}(:,index_fourier{2});
Efrequency(C34((1:Q)+Q)-unknowns,index_space{2}) =  TC4{2}(:,index_fourier{2});

% DC components
Efrequency(dc_index(1)-unknowns,index_space{1})  =  TDCA{1}(index_fourier{1});
Efrequency(dc_index(2)-unknowns,index_space{1})  =  TDCB{1}(index_fourier{1});
Efrequency(dc_index(1)-unknowns,index_space{2})  =  TDCA{2}(index_fourier{2});
Efrequency(dc_index(2)-unknowns,index_space{2})  =  TDCB{2}(index_fourier{2});


%==========================================================================

as1 = diag(h_p{2})*TC1{2}(:,1:end-1) + diag(h_n{2})*TC2{2}(:,1:end-1);
bs1 = diag(h_p{2})*TC3{2}(:,1:end-1) + diag(h_n{2})*TC4{2}(:,1:end-1);

as2 = diag(h_p{2})*TC1{1}(:,1:end-1) + diag(h_n{2})*TC2{1}(:,1:end-1);
bs2 = diag(h_p{2})*TC3{1}(:,1:end-1) + diag(h_n{2})*TC4{1}(:,1:end-1);


ar1 = diag(h_p{1})*TC1{1}(:,1:end-1) + diag(h_n{1})*TC2{1}(:,1:end-1);
br1 = diag(h_p{1})*TC3{1}(:,1:end-1) + diag(h_n{1})*TC4{1}(:,1:end-1);

ar2 = diag(h_p{1})*TC1{2}(:,1:end-1) + diag(h_n{1})*TC2{2}(:,1:end-1);
br2 = diag(h_p{1})*TC3{2}(:,1:end-1) + diag(h_n{1})*TC4{2}(:,1:end-1);



Hs1 = S{2}*as1 + C{2}*bs1 + DC{2}*TDCB{2}(:,1:end-1);
Hs2 = S{2}*as2 + C{2}*bs2 + DC{2}*TDCB{1}(:,1:end-1);
  
Hr1 = S{1}*ar1 + C{1}*br1 + DC{1}*TDCB{1}(:,1:end-1);
Hr2 = S{1}*ar2 + C{1}*br2 + DC{1}*TDCB{2}(:,1:end-1);
  
FourierMatrix = zeros(unknowns);
FourierMatrix(index_space{2},index_space{2}) = ...
                                   Hs1(index_fourier{2},index_fourier{2});

FourierMatrix(index_space{2},index_space{1}) = ...
                                   Hs2(index_fourier{2},index_fourier{1});
                               

FourierMatrix(index_space{1},index_space{1}) = ...
                                   Hr1(index_fourier{1},index_fourier{1});

FourierMatrix(index_space{1},index_space{2}) = ...
                                   Hr2(index_fourier{1},index_fourier{2});
                               
% FourierMatrix1 = zeros(FourierSize);
% 
% FourierMatrix1(1:Q,1:2*Q)                       =  -[diag(p_y1) diag(n_y1)];
% FourierMatrix1((1:Q) + 2*Q,(1:2*Q) + 2*Q)       =  -[diag(p_y1) diag(n_y1)];
% 
% 
% FourierMatrix1(2*Q + (-Q+1:0),[2*Q + (-2*Q+1:0) 4*Q + (-2*Q+1:0)]) = ...
%    -[diag(cos(alpha))*[diag(p_y2) diag(n_y2)], ...
%     -diag(sin(alpha))*[diag(p_y2) diag(n_y2)]];
% 
% FourierMatrix1(4*Q + (-Q+1:0),[2*Q + (-2*Q+1:0) 4*Q + (-2*Q+1:0)]) = ...
%    -[diag(sin(alpha))*[diag(p_y2) diag(n_y2)], ...
%      diag(cos(alpha))*[diag(p_y2) diag(n_y2)]];
% CMatrix = [IA{1} IB{2}*0;IA{1}*0 IB{1}; IA{2}; IB{2}; I];
end