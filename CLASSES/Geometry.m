%==========================================================================
% Geometry.m
% Created: 27.10.2017 - 17:11:51
% By: M. Curti
%
% The Geometry file which takes the geometry info from the XML file and is
% storing it in cells. Also it builds the geometry and metrics
%==========================================================================
classdef Geometry
    %GEOMETRY Builds the geometry according to the geometry xml file
    %   The input of the class is the geometry file
    
    properties
        GeometryElement
        parameters;    % Stored parameters
        points; lines; % Cells with data from points and lines
        elements;
        xi;            % Cells with LGL nodes
        mappings;      % Grid maps from the elements
        metrics;       % Grid metrics from the elements
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = Geometry(filename)
            % Calling the xml file
            A = xmlread(filename);
            obj.GeometryElement = A;
            % Obtaining the parameters
            obj.parameters = xml2matlab(A,'Parameters',0,'Attributes');
            % Loading the parameters in local workspace
            parameter_loader(obj.parameters)
            
            % Evaluating the points
            %--------------------------------------------------------------
            % Loading the number of points
            Np = eval(xml2matlab(A,'Points',0,'p','Attribute'));
            % Evaluating the points formulas and filling the properties of
            % the class
            points = zeros(Np,2);
            for k = 1:Np
                points(k,:) =...
                eval(['[',xml2matlab(A,'point',k-1,...
                                          'address','Attribute'), '];']);
            end
            obj.points.coordinates = points;
            % Evaluating the lines
            %--------------------------------------------------------------
            % Loading the number of lines
            Nl = eval(xml2matlab(A,'Lines',0,'l','Attribute'));
            
            % Types of LGL nodes
            LGLNumber = zeros(Nl,1);
            for k = 1:Nl
                LGLNumber(k) = ...
                       eval(xml2matlab(A,'line',k-1,'degree','Attribute'));
            end
            uniqueLGLNumber = unique(LGLNumber);

            % Building the LGL nodes and quadratures
            xi.nodes  = cell(1,length(uniqueLGLNumber));
            xi.w      = cell(1,length(uniqueLGLNumber));
            xi.D      = cell(1,length(uniqueLGLNumber));

            for k = 1:length(uniqueLGLNumber)
                [csi, w] = ...
                    LegendreGausLobattoNodesAndWeights(uniqueLGLNumber(k));
                 xi.nodes{k} = csi;
                 xi.w{k}     = w;
                 xi.D{k}     = PolynomialDerivativeMatrix(csi);
            end
            xi.line_N   = LGLNumber;
            
            xi.csi_for_all_lines = cell(1,Nl);
            xi.w_for_all_lines   = cell(1,Nl);
            xi.D_for_all_lines   = cell(1,Nl);
            
            for k = 1:Nl
                ii = 1;
                while length(xi.nodes{ii}) ~= (LGLNumber(k)+1)
                      ii = ii + 1;
                end
                xi.csi_for_all_lines{k} = xi.nodes{ii};
                xi.w_for_all_lines{k}   = xi.w{ii};
                xi.D_for_all_lines{k}   = xi.D{ii};
            end
            
            obj.xi = xi;
            % Drawing the lines
            lines = cell(1,Nl);
            grid_lines = cell(1,Nl);
            lines_connectivity  = zeros(Nl,2);
            csi_grid = linspace(-1,1,60);
            for k = 1:Nl
                
                csi = xi.csi_for_all_lines{k};
                % Loading the type of the line
                LineType = xml2matlab(A,'line',k-1,'type','Attribute');

                % Loading the connectivity
                connectivity = eval(['[',xml2matlab(...
                          A,'line',k-1,'connectivity','Attribute'),'];']);
                % Check if the index does not exceed the number of lines
                if max(connectivity)>Np
                    error('Index %d in line %d exceed the number of points %d',...
                        max(connectivity),k,Np)
                end
                      
                      switch LineType
                          case 'segment'
                              [temp_line_x, temp_line_y] = draw_line...
                              (points(connectivity(1),:),...
                               points(connectivity(2),:)...
                                                           ,csi,LineType);
                              [temp_grid_x, temp_grid_y] = draw_line...
                              (points(connectivity(1),:),...
                               points(connectivity(2),:)...
                                                      ,csi_grid',LineType);
                          case 'arcc'
                              % Load the center point
                              Pc = eval(['[',xml2matlab(A,'line',k-1,...
                                  'CenterPoint','Attribute'),'];']);
                              [temp_line_x, temp_line_y] = draw_line...
                              (points(connectivity(1),:),...
                               points(connectivity(2),:),csi,Pc,LineType);
                           
                              [temp_grid_x, temp_grid_y] = draw_line...
                              (points(connectivity(1),:),...
                               points(connectivity(2),:)...
                                                   ,csi_grid',Pc,LineType);
                      end
                lines{k} = [temp_line_x; temp_line_y];
                grid_lines{k} = [temp_grid_x; temp_grid_y];
                lines_connectivity(k,:) = connectivity;
            end
            obj.lines.vector = lines;
            obj.lines.connectivity = lines_connectivity;
            obj.lines.grid = grid_lines;
            
            % Drawing the elements map and calculating the metrics
            % Number of elements
            Nel = eval(xml2matlab(A,'Elements',0,'el','Attribute'));
            mappings.Xm      = cell(1,Nel);
            mappings.Ym      = cell(1,Nel);
            mappings.Xm_grid = cell(1,Nel);
            mappings.Ym_grid = cell(1,Nel);
            metrics.Xcsi     = cell(1,Nel);
            metrics.Ycsi     = cell(1,Nel);
            metrics.Xeta     = cell(1,Nel);
            metrics.Yeta     = cell(1,Nel);
            metrics.J        = cell(1,Nel);
            metrics.W        = cell(1,Nel);
            metrics.w_csi    = cell(1,Nel);
            metrics.w_eta    = cell(1,Nel);
            metrics.Lx       = cell(1,Nel);
            metrics.Ly       = cell(1,Nel);
            metrics.Lpx      = cell(1,Nel);
            metrics.Lpy      = cell(1,Nel);
            metrics.Neu1     = cell(1,Nel);
            metrics.Neu2     = cell(1,Nel);
            metrics.Neu      = cell(2,Nel);
            metrics.N1       = cell(2,Nel);
            metrics.N2       = cell(2,Nel);
            metrics.N        = cell(2,Nel);
            obj.elements.lines  = zeros(Nel,4);
            obj.elements.points = zeros(Nel,4);
            obj.elements.periodic_lines  = zeros(Nel,4);
            obj.elements.periodic_points = zeros(Nel,4);
            obj.elements.periodic = ...
                             eval(['[',xml2matlab(A,'BoundaryConditions'...
                                ,0,'periodic_elements','Attribute'),'];']);
            s      = cell(2,4);
            s_grid = cell(2,4);
            kp = 0; % periodic k
            for k = 1:Nel
                % Loading the connectivities for lines and points
                connectivity = eval(['[',xml2matlab(A,'line_element'...
                                 ,k-1,'connectivity','Attribute'),'];']);
                
                if max(connectivity)>Nl
                    error('Index %d in element %d exceed the number of lines %d',...
                        max(connectivity),k,Nl)
                end
                conn_points = abs(eval(['[',xml2matlab(A,'point_element'...
                                 ,k-1,'connectivity','Attribute'),'];']));
                obj.elements.lines(k,:)  = connectivity; 
                connectivity = abs(connectivity);
                obj.elements.points(k,:) = conn_points;
                if any(k==obj.elements.periodic)
                    obj.elements.periodic_lines(k,:) = ...
                             eval(['[',xml2matlab(A,'periodic_line'...
                                    ,kp,'connectivity','Attribute'),'];']);
                    obj.elements.periodic_points(k,:) = ...
                             eval(['[',xml2matlab(A,'periodic_point'...
                                    ,kp,'connectivity','Attribute'),'];']);
                                kp = kp + 1;
                end
                
                for ii = 1:4
                    s{1,ii}      = lines{connectivity(ii)}(1,:);
                    s{2,ii}      = lines{connectivity(ii)}(2,:);
                    s_grid{1,ii} = grid_lines{connectivity(ii)}(1,:);
                    s_grid{2,ii} = grid_lines{connectivity(ii)}(2,:);
                end
                
                csi   = xi.csi_for_all_lines{connectivity(1)};
                w_csi = xi.w_for_all_lines{connectivity(1)};
                Dcsi  = xi.D_for_all_lines{connectivity(1)};
                
                eta   = xi.csi_for_all_lines{connectivity(2)};
                w_eta = xi.w_for_all_lines{connectivity(2)};
                Deta  = xi.D_for_all_lines{connectivity(2)};
                
                try
                    [temp_Xm, temp_Ym] = TransfiniteQuadMap(csi, eta, s);
                    [temp_Xm_grid, temp_Ym_grid] = ...
                        TransfiniteQuadMap(csi_grid', csi_grid', s_grid);
                    [temp_Xcsi, temp_Ycsi, temp_Xeta, temp_Yeta, temp_J ]...
                        = TransfiniteQuadMetrics(csi, eta, s);
                catch
                    difxi  = numel(s{1,1}) - numel(s{1,3});
                    difeta = numel(s{1,2}) - numel(s{1,4});
                    if difxi > 0
                        error('The lengths of the lines 1 and 3 in the element %d are not consistent',k)
                    elseif  difeta > 0
                        error('The lengths of the lines 2 and 4 in the element %d are not consistent',k)
                    else
                        error('Something is wrong with element %d',k)
                    end
                        
                end

                mappings.Xm{k}      = temp_Xm;
                mappings.Ym{k}      = temp_Ym;
                mappings.Xm_grid{k} = temp_Xm_grid;
                mappings.Ym_grid{k} = temp_Ym_grid;

                metrics.Xcsi{k}  = temp_Xcsi;
                metrics.Ycsi{k}  = temp_Ycsi;
                metrics.Xeta{k}  = temp_Xeta;
                metrics.Yeta{k}  = temp_Yeta;
                metrics.J{k}     = temp_J;
                metrics.W{k}     = w_eta*w_csi';
                metrics.w_csi{k} = w_csi';
                metrics.w_eta{k} = w_eta;
                metrics.Lx{k}    = kron(Dcsi,eye(length(eta)));
                metrics.Ly{k}    = kron(eye(length(csi)),Deta);
                metrics.Lxp{k}   = kron(Dcsi',eye(length(eta)));
                metrics.Lyp{k}   = kron(eye(length(csi)),Deta');
                metrics.Neu1{k}  = sign(temp_J).*sqrt(temp_Yeta.^2 + ...
                                                             temp_Xeta.^2);
                metrics.Neu2{k}  = sign(temp_J).*sqrt(temp_Ycsi.^2 + ...
                                                             temp_Xcsi.^2);
                metrics.N1{1,k}  = sign(temp_J).*temp_Yeta./...
                                         sqrt(temp_Yeta.^2 + temp_Xeta.^2);
                metrics.N1{2,k}  =-sign(temp_J).*temp_Xeta./...
                                         sqrt(temp_Yeta.^2 + temp_Xeta.^2);
                                                       
                
                metrics.N2{1,k} =-sign(temp_J).*temp_Ycsi./...
                                         sqrt(temp_Ycsi.^2 + temp_Xcsi.^2);
                metrics.N2{2,k} = sign(temp_J).*temp_Xcsi./...
                                         sqrt(temp_Ycsi.^2 + temp_Xcsi.^2);
                                     
                % Computing the normal vector on the surface of each
                % element
                % Masks for lines 1,3 and 2,4
                [row, col] = size(temp_Xm); 
                one_row = ones(row,1);
                one_col = ones(1,col);
                
                L13mask = zeros(row, col); L13mask(1,1:end) = -one_col;
                L13mask(end,1:end) = one_col;
                
                
                L24mask = zeros(row, col); L24mask(1:end,1) = one_row;
                L24mask(1:end,end) = -one_row;
                
                metrics.Neu{1,k} = metrics.Neu2{k}.*L13mask*diag(w_csi);
                metrics.Neu{2,k} = diag(w_eta)*L24mask.*metrics.Neu1{k};
                
                metrics.N{1,k} = L13mask.*metrics.N2{1,k} + ...
                                               metrics.N1{1,k}.*(abs(L24mask));
                metrics.N{2,k} = L13mask.*metrics.N2{2,k} + ...
                                               metrics.N1{2,k}.*(abs(L24mask));
                
               tmp = int_el(metrics.N{1,k},'address'); c = tmp{3};
               
               Ncabs = abs(metrics.N{1,k}(c) + metrics.N{2,k}(c)*1i);
               
               metrics.N{1,k}(c) = metrics.N{1,k}(c)./Ncabs;
               metrics.N{2,k}(c) = metrics.N{2,k}(c)./Ncabs;
               
               
            end
            
            % Evaluating the transformations if available
            % Checking the number of transformations
            try
                nT = eval(['[',xml2matlab(A,'Transformations'...
                    ,0,'nT','Attribute'),'];']);
            catch
            end
                               
            if exist('nT','var')
                for t = 1:nT
                    % Extracting the parameters for the transformation
                    type = xml2matlab(A,'transformation'...
                        ,t-1,'type','Attribute');
                    if strcmp(type,'rotation')
                        angle = eval(['[',xml2matlab(A,...
                        'transformation',t-1,...
                        'angle','Attribute'),'];']);
                    end
                    ElementList = eval(['[',xml2matlab(A,...
                        'transformation',t-1,'ElementList','Attribute'),'];']);
                    Copies = eval(['[',xml2matlab(A,...
                        'transformation',t-1,'Copies','Attribute'),'];']);
                    StartLine = eval(['[',xml2matlab(A,...
                        'transformation',t-1,'StartLine','Attribute'),'];']);
                    StartPoint = eval(['[',xml2matlab(A,...
                        'transformation',t-1,'StartPoint','Attribute'),'];']);
                    
                    nEl = numel(ElementList);
                    elements_lines  = obj.elements.lines(ElementList,:);
                    elements_points = obj.elements.points(ElementList,:);
                    
                    inv_lines  = unique(elements_lines(:));
                    inv_points = unique(elements_points(:));
                    
                    num_lines  = numel(inv_lines);
                    num_points = numel(inv_points);
                    % Initialising the new components of the geometry
                    new_el_lines  = zeros(Copies*nEl,4);
                    new_el_points = zeros(Copies*nEl,4);
                    
                    new_mappings.Xm      = cell(1,Copies*nEl);
                    new_mappings.Ym      = cell(1,Copies*nEl);
                    new_mappings.Xm_grid = cell(1,Copies*nEl);
                    new_mappings.Ym_grid = cell(1,Copies*nEl);
                    new_metrics.Xcsi     = cell(1,Copies*nEl);
                    new_metrics.Ycsi     = cell(1,Copies*nEl);
                    new_metrics.Xeta     = cell(1,Copies*nEl);
                    new_metrics.Yeta     = cell(1,Copies*nEl);
                    new_metrics.J        = cell(1,Copies*nEl);
                    new_metrics.W        = cell(1,Copies*nEl);
                    new_metrics.w_csi    = cell(1,Copies*nEl);
                    new_metrics.w_eta    = cell(1,Copies*nEl);
                    new_metrics.Lx       = cell(1,Copies*nEl);
                    new_metrics.Ly       = cell(1,Copies*nEl);
                    new_metrics.Lpx      = cell(1,Copies*nEl);
                    new_metrics.Lpy      = cell(1,Copies*nEl);
                    new_metrics.Neu1     = cell(1,Copies*nEl);
                    new_metrics.Neu2     = cell(1,Copies*nEl);
                    new_metrics.Neu      = cell(2,Copies*nEl);
                    new_metrics.N1       = cell(2,Copies*nEl);
                    new_metrics.N2       = cell(2,Copies*nEl);
                    new_metrics.N        = cell(2,Copies*nEl);
                    for c = 1:Copies
                        new_el_lines((1:nEl) + nEl*(c-1),:) = ...
                            obj.elements.lines + StartLine - 1;
                        StartLine = StartLine + num_lines;
                        
                        new_el_points((1:nEl) + nEl*(c-1),:) = ...
                            obj.elements.points + StartPoint - 1;
                        StartPoint = StartPoint + num_points;
                        
                        for el = (1:nEl) + nEl*(c-1)
                            old_id = el-nEl*(c-1);
                            [TH, R] = cart2pol(mappings.Xm{old_id}...
                                ,mappings.Ym{el-nEl*(c-1)});
                            [Xm, Ym] = pol2cart(TH + angle*pi/180*c, R);
                            
                            new_mappings.Xm{el}   = Xm;
                            new_mappings.Ym{el}   = Ym;
                            new_metrics.Xcsi{el}  = metrics.Xcsi{old_id};
                            new_metrics.Ycsi{el}  = metrics.Ycsi{old_id};
                            new_metrics.Xeta{el}  = metrics.Xeta{old_id};
                            new_metrics.Yeta{el}  = metrics.Yeta{old_id};
                            new_metrics.J{el}     = metrics.J{old_id};
                            new_metrics.W{el}     = metrics.W{old_id};
                            new_metrics.w_csi{el} = metrics.w_csi{old_id};
                            new_metrics.w_eta{el} = metrics.w_eta{old_id};
                            new_metrics.Lx{el}    = metrics.Lx{old_id};
                            new_metrics.Ly{el}    = metrics.Ly{old_id};
                            new_metrics.Lpx{el}   = metrics.Lpx{old_id};
                            new_metrics.Lpy{el}   = metrics.Lpy{old_id};
                            new_metrics.Neu1{el}  = metrics.Neu1{old_id};
                            new_metrics.Neu2{el}  = metrics.Neu2{old_id};
                            new_metrics.Neu{el}   = metrics.Neu{old_id};
                            new_metrics.N1{el}    = metrics.N1{old_id};
                            new_metrics.N2{el}    = metrics.N2{old_id};
                            new_metrics.N{el}     = metrics.N{old_id};
                        end
                    end
                    % Updating the points and lines
                    nE  = Copies*nEl;
                    nP  = max(new_el_points(:));  new_points = zeros(nP,2);
                    nL  = max(abs(new_el_lines(:))); new_lines = cell(1,nL);
                    for el = 1:nE
                        tmp = int_el(new_mappings.Xm{el},'address');
                        p_id = tmp{3}; l_id = tmp{2};
                        
                        xp = new_mappings.Xm{el}(p_id);
                        yp = new_mappings.Ym{el}(p_id);
                        
                        
                        
                        new_points(new_el_points(el,:),:) = [xp' yp'];
                        
                        for ll = 1:4
                            xl = new_mappings.Xm{el}(l_id{ll});
                            yl = new_mappings.Ym{el}(l_id{ll});
                            new_lines{new_el_lines(el,ll)} = [xl; yl];
                        end
                    end
                end
            end
            % Finding the elements that share a line
            shared_line = zeros(Nl,4);
            element_lines = abs(obj.elements.lines);
            for k = 1:Nl
                [rows, cols]= find(k==element_lines);
                for r = 1:length(rows)
                    shared_line(k,r)   = rows(r);
                    shared_line(k,r+2) = ...
                        cols(r)*sign(obj.elements.lines(rows(r),cols(r)));
                end
            end
            
            % Finding the elements that share a point
            shared_point = cell(2,Nl);
            for k = 1:Np
                [rows, cols] = find(k == obj.elements.points);
                
                shared_point{1,k} = rows;
                shared_point{2,k} = cols;
            end
            
            obj.lines.shared_line   = shared_line;
            obj.points.shared_point = shared_point;
            
            obj.mappings = mappings;
            obj.metrics  = metrics;
            
        end
        %==================================================================

        %------------------------------------------------------------------
        % Function to plot geometry entities
        %------------------------------------------------------------------
        function plot_geometry(obj,varargin)
            mode = varargin{nargin-1};
            switch mode
               % Plot points
               case 'points'
               PointsNumber = eval(xml2matlab(obj.GeometryElement...
                                            ,'Points',0,'p','Attribute'));
               if nargin > 2
                  PropertiesNumber = (nargin-2)/2;
               end
               for k = 1:PointsNumber
                   hp = plot(obj.points.coordinates(k,1),...
                             obj.points.coordinates(k,2),'.');
                         text(obj.points.coordinates(k,1)-1,obj.points.coordinates(k,2),sprintf('%d',k),'HorizontalAlignment','right','Color','red')

                   if exist('PropertiesNumber','var')
                      for ii = 1:PropertiesNumber
                          set(hp,varargin{ii*2-1},varargin{ii*2})
                      end
                   end
               end
               % Plot lines
               
               case 'lines'
               LinesNumber = eval(xml2matlab(obj.GeometryElement...
                                            ,'Lines',0,'l','Attribute'));
               if nargin > 2
                  PropertiesNumber = (nargin-2)/2;
               end
               for k = 1:LinesNumber
                   temp_line = obj.lines.vector{k};
                   hp = plot(temp_line(1,:),temp_line(2,:));
                   middle = round(numel(temp_line(1,:))/2);
%                    text(temp_line(1,middle)-1,temp_line(2,middle)-1,sprintf('%d',k),'HorizontalAlignment','right','Color','green')

                   if exist('PropertiesNumber','var')
                      for ii = 1:PropertiesNumber
                          set(hp,varargin{ii*2-1},varargin{ii*2})
                      end
                   end
               end
               % Plot element grid
               
               case 'ElementGrid'
               ElementNumber = ...
                        eval(xml2matlab(obj.GeometryElement,...
                                          'Elements',0,'el','Attribute'));
               if nargin > 2
                  PropertiesNumber = (nargin-2)/2;
               end
               for k = 1:ElementNumber
                   hp = plot(obj.mappings.Xm{k},obj.mappings.Ym{k},'.');

                   if exist('PropertiesNumber','var')
                      for ii = 1:PropertiesNumber
                          set(hp,varargin{ii*2-1},varargin{ii*2})
                      end
                   end
               end

            end
        end
        
        
        %------------------------------------------------------------------
        % Update Parameters
        %------------------------------------------------------------------
        function obj = setParameter(obj,Parameter,Value)
            obj = change_parameter( obj,Parameter,Value );
        end
    end
    
end