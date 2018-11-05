function [ Ic0, Ia, Ib ] = lagrange2fourier_int( x,Q )
%LAG2FOURIER_INT returns the matrices which, multiplied by the functions on
%lag2fourier Transfomation into sine and cosine coefficients from
%polynomial nodes.
%   xi_j - abscisas of the polynomial nodes
%   w_j  - quadratures or any operator weights
%   H    - number of harmonics

for ii = 1:50
    Qk = ii*6;
    if not(iscell(x))
        
        %==================================================================
        % Routine for a single element
        %------------------------------------------------------------------
        % Reshaping the vectors
        x = reshape(x,1,numel(x));
        %------------------------------------------------------------------
        
        tau   = (x(end) - x(1))*.5; % The period of the function
        w_i   = (1:Q)'*pi/tau;      % Calculated frequencies based on Q
        
        % Preparing the interpolation basis functions and the quadratures
        
        M     = Qk;                % Size of the integration vector
        [x_int, w_int] = LegendreGausLobattoNodesAndWeights(M);
        
        % Scalling the integration nodes and weights to the limits of x
        x_int = n2range(x_int,x(1),x(end));
        w_int = w_int*(x(end)-x(1))/2;
        
        % Calculating the interpolation basis function at integration nodes
        l = PolynomialInterp(x_int,x);
        
        % Expanding into Fourier series
        % DC component
        Ic0   = ((l'*w_int)/(2*tau))';
        
        % sine components
        Ia    = sin(w_i*x_int')*diag(w_int)*l/tau;
        
        % cosine components
        Ib    = cos(w_i*x_int')*diag(w_int)*l/tau;
        
        
    else
        % Routine for multiple elements stored in a cell
        %==================================================================
        
        % Extracting the problem data
        Nel       = numel(x);     % Number of elements
        Nnodes    = zeros(1,Nel); % Number of nodes in each element
        tau       = (x{end}(end) - x{1}(1))/2;
        % Period of the function
        w_i       = (1:Q)'*pi/tau;% Calculated frequencies based on Q
        
        
        % Evaluate the basis functions per interval
        
        Ic0c = cell(1,Nel); Iac = cell(1,Nel); Ibc = cell(1,Nel);
        
        for k = 1:Nel
            Nnodes(k) = numel(x{k});
            M         = Nnodes(k)*ii;          % Size of the integration vector
            % Calculating the nodes and weights for the integration
            [x_int, w_int] = LegendreGausLobattoNodesAndWeights(M);
            
            % Scaling the nodes and weights
            x_int_n = n2range(x_int,x{k}(1),x{k}(end));
            w_int_n = w_int*(x{k}(end)-x{k}(1))/2;
            
            % Evaluate the interpolation basis at the integration nodes
            l       = PolynomialInterp(x_int_n,x{k});
            
            % Evaluate the Fourier coefficients for each element
            % The DC component
            Ic0c{k} = ((l'*w_int_n)/(2*tau))';%reshape(w{k},1,numel(w{k}))/(2*tau);%
            
            % sine components
            Iac{k}  = sin(w_i*x_int_n')*diag(w_int_n)*l/tau;
            
            % cosine components
            Ibc{k}  = cos(w_i*x_int_n')*diag(w_int_n)*l/tau;
        end
        NnodesSum = sum(Nnodes);  % Total number of nodes;
        
        
        % Grouping the cell integrands into one matrix
        Ic0 = zeros(1,NnodesSum - Nel + 1); 
        Ia  = zeros(Q,NnodesSum - Nel + 1);
        Ib  = zeros(Q,NnodesSum - Nel + 1);
        
        currentIndex = 0;
        for k = 1:Nel
            if k == 1
                Ic0(1:Nnodes(k))  = Ic0c{k};
                Ia(:,1:Nnodes(k)) = Iac{k};
                Ib(:,1:Nnodes(k)) = Ibc{k};
                
                currentIndex = Nnodes(k);
            else
                Ic0((0:Nnodes(k)-1) + currentIndex) = ...
                          Ic0((0:Nnodes(k)-1) + currentIndex) + Ic0c{k};
                Ia(:,(0:Nnodes(k)-1) + currentIndex)  = ...
                          Ia(:,(0:Nnodes(k)-1) + currentIndex) + Iac{k};
                Ib(:,(0:Nnodes(k)-1) + currentIndex)  = ...
                          Ib(:,(0:Nnodes(k)-1) + currentIndex) + Ibc{k};
                
                
                currentIndex = currentIndex + Nnodes(k) - 1;
            end
        end
        
    end
    
    % Evaluating the error of the interpolands
    
    if ii == 1
        prevIc0 = Ic0; prevIa = Ia; prevIb = Ib;
    else
        dIc0 = prevIc0 - Ic0; dIa = prevIa - Ia; dIb = prevIb - Ib;
        
        max_err = max(abs([dIc0(:); dIa(:); dIb(:)]));
        prevIc0 = Ic0; prevIa = Ia; prevIb = Ib;
    end
    if ii > 1 && max_err < 3e-15
        break
    end
end

end
