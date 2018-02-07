%==========================================================================
% fourierElements.m
% Created: 03.11.2017 - 10:11:59
% By: M. Curti
%
%fourierElements Performs the computations on one domain and returns matrix
%components
%==========================================================================
classdef fourierElements
    %FOURIER Performs the computations on one domain and returns matrix
    %components
    %   Detailed explanation goes here
    
    properties
        Edata   % Data of the Fourier Elements
        
        ys   % Scaling factors for 
        w_n  % The frequencies for each harmonic
        interface_line
        s; c1; c2; c3; c4; Az0; Bx0;
             % Coefficients for the solution
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        function obj = fourierElements(Elements)
            
            obj.Edata      = Elements;
            
            obj.w_n        = 2*(1:Elements.Harmonics)'*pi/Elements.tau;
            obj.Edata.Nfel = numel(Elements.heights)/2;
            
            % Computing the scalling constants
            
            ys = zeros(obj.Edata.Nfel,2);
            
            for k = 1:obj.Edata.Nfel
                ys(k,1) = max(Elements.heights(:));
                ys(k,2) = min(Elements.heights(:));
            end
            obj.ys = ys;
        end
        %------------------------------------------------------------------
        
        % Submatrix for the Vector magnetic potential
        function [azs, azc] = Az_builder(obj, y)
                        
            azs = [diag(exp(obj.w_n*(y-obj.ys1)))...
                                          diag(exp(-obj.w_n*(y-obj.ys2)))];
            azc = [diag(exp(obj.w_n*(y-obj.ys1)))...
                                          diag(exp(-obj.w_n*(y-obj.ys2)))];
        end
        
        %------------------------------------------------------------------
        % Function to update the solutions coefficients
        function obj = update_coefficients(obj, s, c1, c2, c3, c4, Az0, Bx0)
            obj.c1 = c1; obj.c2 = c2; obj.c3 = c3; obj.c4 = c4; obj.s = s;
            obj.Az0 = Az0; obj.Bx0 = Bx0;
        end
        
        %------------------------------------------------------------------
        % Function to evaluate the spatial values of the solution within
        % the defined domain
        function [x_grid, y_grid, f_solution] = ...
                                             fourier2space(obj, dx, dy, El)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
            
            % Intialisation
            Q  = obj.Edata.Harmonics;
            x1 = obj.Edata.x_start;       x2 = x1 + obj.Edata.tau;
            y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);
            C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
            C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
            C0 = obj.s((1:2*Q) + (El-1)*2*Q);
            %==============================================================
            
            x_i = linspace(x1, x2, dx); y_j = linspace(y1, y2, dy);
            
            [x_grid, y_grid] = meshgrid(x_i,y_j);
            
            a  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1) + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2) - repmat(C0((1:Q)+Q)./obj.w_n',dy,1);
            b  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3) + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4) + repmat(C0(1:Q)./obj.w_n',dy,1);
             
            c0 = obj.Az0 + obj.Bx0*y_grid;
            f_solution = c0 + a*sin(obj.w_n*x_i) + b*cos(obj.w_n*x_i);
        end
        
        %------------------------------------------------------------------
        % Function to evaluate the spatial values for the derivative of the
        % solution with the respect to y axis, namely Bx
        function [x_grid, y_grid, f_solution] = ...
                fourier_Bx2space(obj, dx, dy,El)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
            
            % Intialisation
            Q  = obj.Edata.Harmonics;
            x1 = obj.Edata.x_start;       x2 = x1 + obj.Edata.tau;
            y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);
            C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
            C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
%             C0 = obj.s((1:Q) + (El-1)*Q);
            %==============================================================
            x_i = linspace(x1, x2, dx);
            y_j = linspace(y1, y2, dy);
            
            [x_grid, y_grid] = meshgrid(x_i,y_j);
            
            a  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') - ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n');
            b  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') - ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n');
             
            c0 = obj.Bx0;
            f_solution = c0 + a*sin(obj.w_n*x_i) + b*cos(obj.w_n*x_i);
            
        end
        
        %------------------------------------------------------------------
        % Function to evaluate the spatial values for the derivative of the
        % solution with the respect to y axis, namely Bx
        function [x_grid, y_grid, f_solution] = ...
                                          fourier_By2space(obj, dx, dy, El)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
            % Intialisation
            Q  = obj.Edata.Harmonics;
            x1 = obj.Edata.x_start;       x2 = x1 + obj.Edata.tau;
            y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);
            C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
            C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
            C0 = obj.s((1:2*Q) + (El-1)*2*Q);
            %==============================================================
            x_i = linspace(x1, x2, dx);
            y_j = linspace(y1, y2, dy);
            
            [x_grid, y_grid] = meshgrid(x_i,y_j);
            
            a  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n') + repmat(C0((1:Q)+Q).*obj.w_n',dy,1) ;
            b  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n') - repmat(C0(1:Q).*obj.w_n',dy,1) ;
             
            c0 = 0;
            f_solution = -c0 - a*cos(obj.w_n*x_i) + b*sin(obj.w_n*x_i);
        end
        %------------------------------------------------------------------
        % Function to plot Fourier solution on mesh nodes
        function f_solution = fourier2mesh(obj, x, y,El)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
             % Intialisation
            Q  = obj.Edata.Harmonics;
            C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
            C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
            C0 = obj.s((1:2*Q) + (El-1)*2*Q);
            %==============================================================
            points_number = numel(x);
            
            f_solution = zeros(points_number,1);
            
            for k = 1:points_number
                a  = exp( obj.w_n*(y(k)-obj.ys(El,1))).*C1' + ...
                     exp(-obj.w_n*(y(k)-obj.ys(El,2))).*C2'+ C0((1:Q)+Q)'./obj.w_n;
                b  = exp( obj.w_n*(y(k)-obj.ys(El,1))).*C3' + ...
                     exp(-obj.w_n*(y(k)-obj.ys(El,2))).*C4'+ C0(1:Q)'./obj.w_n;
                c0 = obj.Az0 + obj.Bx0*y(k);
                
                f_solution(k) = c0 + sum(a.*sin(obj.w_n*x(k)) + ...
                                                     b.*cos(obj.w_n*x(k)));
            end
        end
        
        %------------------------------------------------------------------
        % Function to plot Fourier solution on mesh nodes
        function f_solution = fourier_Bx2mesh(obj, x, y)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
            
            points_number = numel(x);
            
            f_solution = zeros(points_number,1);
            
            for k = 1:points_number
                a  = exp( obj.w_n*(y(k)-obj.ys1)).*obj.c1.*obj.w_n - ...
                             exp(-obj.w_n*(y(k)-obj.ys2)).*obj.c2.*obj.w_n;
                b  = exp( obj.w_n*(y(k)-obj.ys1)).*obj.c3.*obj.w_n - ...
                             exp(-obj.w_n*(y(k)-obj.ys2)).*obj.c4.*obj.w_n;
                c0 = obj.Bx0;
                
                f_solution(k) = c0 + sum(a.*sin(obj.w_n*x(k)) + ...
                                                     b.*cos(obj.w_n*x(k)));
            end
        end
        %------------------------------------------------------------------
        % Function to plot Fourier solution on mesh nodes
        function f_solution = fourier_By2mesh(obj, x, y)
            
            if isempty(obj.c1)
                error('The coefficients for the solution are not upladed!')
            end
            
            points_number = numel(x);
            
            f_solution = zeros(points_number,1);
            
            for k = 1:points_number
                a  = exp( obj.w_n*(y(k)-obj.ys1)).*obj.c1 + ...
                                      exp(-obj.w_n*(y(k)-obj.ys2)).*obj.c2;
                b  = exp( obj.w_n*(y(k)-obj.ys1)).*obj.c3 + ...
                                      exp(-obj.w_n*(y(k)-obj.ys2)).*obj.c4;
                c0 = 0;
                
                f_solution(k) = -c0 - sum(a.*cos(obj.w_n*x(k)).*obj.w_n-...
                                            b.*sin(obj.w_n*x(k)).*obj.w_n);
            end
        end
        %------------------------------------------------------------------
        % Function to build the part of the Global matrix which coresponds
        % to Fourier evaluated in spatial domain. For this function the
        % ProblemData structure is needed at the input
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
            for k = 1:Nel
               Nnodes(k)     = ...
                          numel(Geometry.lines.vector{conn_lines(k)}(1,:));
            end
            NnodesSum = sum(Nnodes);  % Total number of nodes;
            
            % Building the quadrature and abscisas vector
            
            x_fourier     = zeros(1,NnodesSum-Nel+1);
            w_fourier     = zeros(1,NnodesSum-Nel+1);
            
            currentIndex = Nnodes(1);
            for k = 1:Nel
                %------------------------------------------------------
                neul       = int_el(metrics.Neu2{conn_El(k)},1,'lines');
                neup       = int_el(metrics.Neu2{conn_El(k)},...
                                                           [1 2],'points');
                neu2       = [neup(1) neul{1} neup(2)];
                l = conn_lines(k);
                if k==1
                    x_fourier(1:Nnodes(1)) = Geometry.lines.vector{l}(1,:);
                    w_fourier(1:Nnodes(1)) = ...
                                Geometry.xi.w_for_all_lines{l}'.*neu2;                     
                else
                    x_fourier((0:Nnodes(k)-1) + currentIndex) = ...
                                             Geometry.lines.vector{l}(1,:);
                    w_fourier((0:Nnodes(k)-1) + currentIndex) = ...
                            w_fourier((0:Nnodes(k)-1) + currentIndex) + ...   
                                Geometry.xi.w_for_all_lines{l}'.*neu2; 
                    currentIndex = currentIndex + Nnodes(k) - 1;
                end
            end
            % Removing the last entry in the vectors
%             x_fourier(end) = []; w_fourier(end) = [];
            
            % Initiating the Espace matrix
            unknowns = ProblemData.inElementsUnknowns + ...
                       ProblemData.linesUnknowns + ...
                       ProblemData.pointsUnknowns;
            Espace   = zeros(unknowns,Q*4*Nfel + 1);
                                    
            % Initiating the entries, the exponentials are necessary only
            % for the first layer, since other layers are not connected to
            % the SEM part
            
            [p_exp, n_exp] = np_exp(obj.w_n,obj.Edata.heights(1,1),...
                                                  obj.ys(1,1),obj.ys(1,2));
                                              
            w_fourier([1 end]) = w_fourier(1) + w_fourier(end);

            % Computing the block matrix in the Espace matrix
            Espace_block = -diag(w_fourier(1:end-1))*...
                 [  sin(x_fourier(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
                  - sin(x_fourier(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
                    cos(x_fourier(1:end-1)'*obj.w_n')*diag(p_exp)*diag(obj.w_n), ...
                  - cos(x_fourier(1:end-1)'*obj.w_n')*diag(n_exp)*diag(obj.w_n), ...
                    zeros(numel(x_fourier(1:end-1)),Q*4*(Nfel-1)),...
                     ones(numel(x_fourier(1:end-1)),1)];
                 
%             Espace_block = -diag(w_fourier)*...
%                  [- cos(x_fourier'*obj.w_n')*diag(p_exp), ...
%                     cos(x_fourier'*obj.w_n')*diag(n_exp), ...
%                     sin(x_fourier'*obj.w_n')*diag(p_exp), ...
%                   - sin(x_fourier'*obj.w_n')*diag(n_exp), ...
%                     zeros(numel(x_fourier),Q*4*(Nfel-1)),...
%                      x_fourier'];

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
                    f_lines((1:Nnodes(k)-2) + lineIndex) = ...
                                   tmp_inc((2:Nnodes(k)-1) + currentIndex);
                    currentIndex = currentIndex + Nnodes(k)-1;
                    lineIndex    = lineIndex + Nnodes(k)-2;
                end
            end
            index_fourier = [f_lines P(1:end-1)];
            % Moving the block to the matrix
            Espace(index_space,:) = Espace_block(index_fourier,:);
            
            
            interface_line.x = x_fourier;
            interface_line.w = w_fourier;
        end
        %------------------------------------------------------------------
        % Function to build the part of the Global matrix which coresponds
        % to Fourier sapce. For this function the
        % ProblemData structure is needed at the input
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
           
           %===============================================================
           for k = 1:Nel
               x_fourier{k}  = Geometry.lines.vector{conn_lines(k)}(1,:);
               Nnodes(k)     = numel(x_fourier{k});
               l = conn_El(k);
               %------------------------------------------------------
               neul       = int_el(metrics.Neu2{conn_El(l)},1,'lines');
               neup       = int_el(metrics.Neu2{conn_El(l)},...
                                                           [1 2],'points');
               %------------------------------------------------------
               neu2       = [neup(1) neul{1} neup(2)];
               w_fourier{k}  = ...
                         Geometry.xi.w_for_all_lines{conn_lines(k)}'.*neu2;
           end
           NnodesSum = sum(Nnodes);  % Total number of nodes;
           % SEM to fourier transformation
           [Ic0, Ia, Ib]   = lagrange2fourier_int(x_fourier,w_fourier, Q);
%            [Ic0, Ia, Ib]   = lagrange2fourier(x_fourier,w_fourier,Q);
            Ic0(1) = Ic0(1)+Ic0(end);
            Ia(:,1) = Ia(:,1)+Ia(:,end);
            Ib(:,1) = Ib(:,1)+Ib(:,end);
           % Initiating the Efrequency matrix
            unknowns       = ProblemData.inElementsUnknowns + ...
                             ProblemData.linesUnknowns + ...
                             ProblemData.pointsUnknowns; % The unknowns in 
                                                         % the SEM model
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
                    f_lines((1:Nnodes(k)-2) + lineIndex) = ...
                                   tmp_inc((2:Nnodes(k)-1) + currentIndex);
                    currentIndex = currentIndex + Nnodes(k)-1;
                    lineIndex    = lineIndex + Nnodes(k)-2;
                end
            end
           index_fourier = [f_lines P(1:end-1)];


           % Building the boundary condition matrix
           
           FourierMatrix = zeros(FourierSize);
           K = @(k) k*4*Q;
           for k = 1:Nfel
               [p_exp_y1, n_exp_y1] = ...
                           np_exp(obj.w_n,hf(k,1),obj.ys(k,1),obj.ys(k,2));
               [p_exp_y2, n_exp_y2] = ...
                           np_exp(obj.w_n,hf(k,2),obj.ys(k,1),obj.ys(k,2));
                              
               if k == 1
                   FourierMatrix(1:Q,1:2*Q) = ...
                           -[diag(p_exp_y1)  diag(n_exp_y1)];
                   FourierMatrix((1:Q) + 2*Q,(1:2*Q) + 2*Q) = ...
                           -[diag(p_exp_y1)  diag(n_exp_y1)];
                   % The DC component with SEM
%                    FourierMatrix(Nfel*Q*4 + 1,Nfel*Q*4 + 1) = ...
%                                                             -hf(k,1);
               elseif k == 2
                   %Az
                   %sine
                   FourierMatrix((1:Q)+Q,(1:6*Q)) = ...
                     [diag(p_exp_y1),   diag(n_exp_y1), ...
                       zeros(Q,2*Q), -diag(p_exp_y1), -diag(n_exp_y1)];
                   %cosine
                   FourierMatrix((1:Q)+3*Q,(1:6*Q)+2*Q) = ...
                     [diag(p_exp_y1) diag(n_exp_y1), ...
                       zeros(Q,2*Q), -diag(p_exp_y1), -diag(n_exp_y1)];
                   
                   %Hx
                   %sine
                   FourierMatrix((1:Q)+4*Q,(1:6*Q)) = ...
                       [diag(p_exp_y1), -diag(n_exp_y1), ...
                       zeros(Q,2*Q), -diag(p_exp_y1),   diag(n_exp_y1)];
                   
                   %cosine
                   FourierMatrix((1:Q)+6*Q,(1:6*Q)+2*Q) = ...
                       [diag(p_exp_y1), -diag(n_exp_y1), ...
                       zeros(Q,2*Q), -diag(p_exp_y1),   diag(n_exp_y1)];
                   % The DC components to space
%                    FourierMatrix(end-1,[K(k)+2 end]) =  [hf(k,1) -hf(k,1)];
               end
               
               
               
               if k == Nfel
                   FourierMatrix(K(k) - 2*Q + (-Q+1:0), K(k) - ...
                       2*Q + (-2*Q+1:0)) = ...
                            [diag(p_exp_y2) -diag(n_exp_y2)];
                         
                   FourierMatrix(K(k) + (-Q+1:0),K(k) + (-2*Q+1:0)) = ...
                            [diag(p_exp_y2) -diag(n_exp_y2)];
                   % The DC components to space
                   FourierMatrix(end,K(k) + 1) = -hf(k,1);
               end
               
           end
           % Building the Efrequency based on boundary conditions

           % Sine term for the shared line
           % SEM side
           Efrequency(sin_index(1:end/2)-unknowns,index_space) = ...
                                         Ia(:,index_fourier);
           % Cosine terms
           % SEM side
           Efrequency(cos_index(1:end/2)-unknowns,index_space) = ...
                                        Ib(:,index_fourier);
           % DC components
           Efrequency(dc_index(1)-unknowns,index_space)        = ...
                                         Ic0(index_fourier);
           Efrequency = [Efrequency FourierMatrix];
           
        end
        
    %------------------------------------------------------------------
    % Function to evaluate the forces
        function Fx = MaxwellForce(obj,El)
           if isempty(obj.c1)
               error('The coefficients for the solution are not upladed!')
           end
           
           % Intialisation
           Q  = obj.Edata.Harmonics;
           
           y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);
           C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
           C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
           C0 = obj.s((1:Q) + (El-1)*Q);
           %==============================================================           x_i = linspace(x1, x2, dx);
           y_j = y1 + (y2-y1)/2;
           
                      
           ax  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') - ...
               exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n');
           bx  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') - ...
               exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n');
           
           
           ay  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C1.*obj.w_n') + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C2.*obj.w_n');
           by  = exp( obj.w_n*(y_j-obj.ys(El,1)))'*diag(C3.*obj.w_n') + ...
                 exp(-obj.w_n*(y_j-obj.ys(El,2)))'*diag(C4.*obj.w_n') - ...
                 C0.*obj.w_n';
           
           Fx = (sum(ax.*by) -  sum(bx.*ay))/(pi*4e-7)*100/obj.w_n(1)*pi;
           
        end
    %------------------------------------------------------------------
    % Function to update the x_start. Used for moving regions
        function obj = update_x_start(obj, new_x)
            obj.Edata.x_start = new_x;
        end
    %------------------------------------------------------------------
    % Function to compute the force in a Fourier region based on wirtual
    % wor
        function Fx = virtual_force(obj, x1, x2, El)
           % Intialisation
           Q  = obj.Edata.Harmonics;
           wn = obj.w_n';
           y1 = obj.Edata.heights(El,2); y2 = obj.Edata.heights(El,1);
           C1 = obj.c1((1:Q) + (El-1)*Q);   C2 = obj.c2((1:Q) + (El-1)*Q);
           C3 = obj.c3((1:Q) + (El-1)*Q);   C4 = obj.c4((1:Q) + (El-1)*Q);
           m  = obj.s((1:Q) + (El-1)*Q);
           mu0 = pi*4e-7;
           % Initial precomputations
           ep = exp(2*wn*(y1-obj.ys(El,1)))-exp(2*wn*(y2-obj.ys(El,1)));
           em = exp(-2*wn*(y1-obj.ys(El,2)))-exp(-2*wn*(y2-obj.ys(El,2)));
           
           ep1 = exp(wn*(y1-obj.ys(El,1)))-exp(wn*(y2-obj.ys(El,1)));
           em1 = exp(-wn*(y1-obj.ys(El,2)))-exp(-wn*(y2-obj.ys(El,2)));
           
           Ixya = 2*C1.^2.*ep + 2*C2.^2.*em - 4*m.*(C1.*ep1 + C2.*em1) +...
                  2*wn*(y2-y1).*(m.^2 + C1.*C2);
           Ixyb = 2*C3.^2.*ep + 2*C4.^2.*em + 4*m.*(C3.*ep1 + C4.*em1) +...
                  2*wn*(y2-y1).*(m.^2 + C3.*C4);
           Ixyc = 2*C1.*C2.*ep + 2*C2.*C4.*em + 2*(C1 - C3).*ep1.*m + ...
                  2*(C2 - C4).*em1.*m + 2*wn*(y2-y1).*m.^2;
           
           % Computation of the forces
           
           Fx = 100/(2*mu0)*sum(wn.*(...
                              (sin(wn.*x2).^2 - sin(wn.*x1).^2).*Ixya + ...
                              (cos(wn.*x2).^2 - cos(wn.*x1).^2).*Ixyb + ...
                              (sin(2*wn.*x2) - cos(2*wn.*x1)).*Ixyc*.5 ...
               ));
        end
    
    end
    
end


% the set of linearly independent functions
function [p_exp, n_exp] = np_exp(w_n,y,s1,s2)
    p_exp = exp( w_n*(y-s1));
    n_exp = exp(-w_n*(y-s2));
end