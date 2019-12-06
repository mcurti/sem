%==========================================================================
% fourierElements.m
% Created: 03.11.2017 - 10:11:59
% By: M. Curti
%
%fourierElements Performs the computations on one domain and returns matrix
%components
%==========================================================================
classdef FourierElement
    %FOURIER Performs the computations on one domain and returns matrix
    %components
    %   Detailed explanation goes here
    
    properties
        Edata   % Data of the Fourier Element
        SEMdata % Data of the SEM elements
        ys   % Scaling factors for
        w_n  % The frequencies for each harmonic
        interface_line
        s; c1; c2; c3; c4; Az0; Bx0;
        % Coefficients for the solution
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        function obj = FourierElement(xmlContent,Geometry)
            % Geting general parameters for Fourier
            Q     = eval(xml2matlab(xmlContent...
                ,'Fourier',0,'Harmonics','Attribute'));
            type  = xml2matlab(xmlContent...
                ,'Fourier',0,'type','Attribute');
            start = eval(xml2matlab(xmlContent...
                ,'Fourier',0,'start','Attribute'));
            tau   = eval(xml2matlab(xmlContent...
                ,'Fourier',0,'tau','Attribute'));
            %                 fel   = eval(xml2matlab(xmlContent...
            %                     ,'Fourier',0,'fel','Attribute'));
            
            % Geting parameters per domain
            parameters = xml2matlab(xmlContent,'Parameters',0,'Attributes');
            parameter_loader(parameters)
            domain.heights    = zeros(1,2);
            domain.heights(1,:) = eval(['[' xml2matlab(xmlContent...
                ,'domain',0,'heights','Attribute') '];']);
            domain.lines_top = eval(['[' xml2matlab(xmlContent...
                ,'domain',0,'lines_top','Attribute') '];']);
            domain.elements_top = eval(['[' xml2matlab(xmlContent...
                ,'domain',0,'elements_top','Attribute') '];']);
            
            domain.lines_bot = eval(['[' xml2matlab(xmlContent...
                ,'domain',0,'lines_bot','Attribute') '];']);
            domain.elements_bot = eval(['[' xml2matlab(xmlContent...
                ,'domain',0,'elements_bot','Attribute') '];']);
            
            % Formating the data for Fourier class
            Elements           = struct;
            Elements.x_start   = start;
            Elements.tau       = tau;
            Elements.Harmonics = Q;
            Elements.heights   = domain.heights;
            Elements.type      = type;
            Geom_data.lines    = Geometry.lines;
            Geom_data.xi       = Geometry.xi;
            Geom_data.metrics  = Geometry.metrics;
            %
            
            Geom_data.Lines_top    = domain.lines_top;
            Geom_data.Elements_top = domain.elements_top;
            Geom_data.Lines_bot    = domain.lines_bot;
            Geom_data.Elements_bot = domain.elements_bot;
            
            obj.Edata      = Elements;
            obj.SEMdata    = Geom_data;
            obj.w_n        = 2*(1:Elements.Harmonics)'*pi/Elements.tau;
            obj.Edata.Nfel = numel(Elements.heights)/2;
            try
                obj.Edata.type = Elements.type;
            catch
                obj.Edata.type = 'cartesian';
                disp('Default coordinate system - cartesian is selected');
            end
            
            % Computing the scalling constants
            
            y_s = zeros(1,2);
            
            for k = 1:obj.Edata.Nfel
                y_s(k,1) = max(Elements.heights(:));
                y_s(k,2) = min(Elements.heights(:));
            end
            obj.ys = y_s;
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
        [x_grid, y_grid, f_solution] = fourier2space(obj, dx, dy, El)
        
        
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
        
        f_solution = fourier2mesh(obj, x, y,El)
        
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
        [Espace, interface_line] = fourier_space_matrix(obj,...
            ProblemData,Geometry,FourierData)
        %------------------------------------------------------------------
        % Function to build the part of the Global matrix which coresponds
        % to Fourier sapce. For this function the
        % ProblemData structure is needed at the input
        Efrequency = fourier_frequency_matrix...
            (obj,ProblemData,Geometry,FourierData)
        % For doffree
        [Efrequency, CM] = fourier_frequency_matrix_doffree...
            (obj,ProblemData,Geometry,FourierData)
        %------------------------------------------------------------------
        % Function to evaluate the forces
        Fx = MaxwellForce(obj,El)
        
        %------------------------------------------------------------------
        % Function to update the x_start. Used for moving regions
        function obj = update_x_start(obj, new_x)
            obj.Edata.x_start = new_x;
        end
        
        %------------------------------------------------------------------
        % Function to develop the doffree method
        tmp = doffree(obj)
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