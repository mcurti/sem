%==========================================================================
% Problem.m
% Created: 27.10.2017 - 14:45:01
% By: M. Curti
%
% Building the global matrix and the right hand side source vector
%==========================================================================
classdef Problem
    %PROBLEM Builds the matrices for the problem
    %   Detailed explanation goes here
    
    properties
        Y;             % Right hand side of the equation
        ProblemData;   % Various parameters of the problem
        Sources;       % Sources in the geometry
        metrics        % Metrics of the geometry
        xmlContent     % Info from the xml file
        Materials      % Materials of the problem
        tmp_Element_Matrix
        Global_Matrix
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = Problem(dataProblem)
            obj.ProblemData = dataProblem.ProblemData;
            obj.Sources     = dataProblem.Sources;
            obj.metrics     = dataProblem.metrics;
            obj.xmlContent  = dataProblem.xmlContent;
            obj.Materials   = dataProblem.Materials;
        end
        %------------------------------------------------------------------
        % Load the Y sources, right hand side
        %------------------------------------------------------------------
        
        function obj = load_Y_sources(obj)
            Nel    = obj.ProblemData.Nel;
            
            obj.Y.CurrentDensity = cell(1,Nel);
            obj.Y.Magnetisation  = cell(1,Nel);
                        
            for k = 1:Nel
                obj.Y.CurrentDensity{k} = -obj.Sources.CurrentDensity{k}...
                    .*obj.metrics.J{k}.*obj.metrics.W{k};
                
                %==========================================================
                
                dataElement.Lx    = obj.metrics.Lx{k};
                dataElement.Ly    = obj.metrics.Ly{k};
                dataElement.Lxp   = obj.metrics.Lxp{k};
                dataElement.Lyp   = obj.metrics.Lyp{k};
                dataElement.Xcsi  = obj.metrics.Xcsi{k};
                dataElement.Ycsi  = obj.metrics.Ycsi{k};
                dataElement.Xeta  = obj.metrics.Xeta{k};
                dataElement.Yeta  = obj.metrics.Yeta{k};
                dataElement.Neu1  = obj.metrics.Neu1{k};
                dataElement.Neu2  = obj.metrics.Neu2{k};
                dataElement.J     = obj.metrics.J{k};
                dataElement.w_csi = obj.metrics.w_csi{k};
                dataElement.w_eta = obj.metrics.w_eta{k};
                dataElement.nu    = 1./obj.Materials.Permeability{k};
                dataElement.Mx    = obj.Sources.Magnetisation{1,k};
                dataElement.My    = obj.Sources.Magnetisation{2,k};
                
                obj.Y.Magnetisation{1,k} = magnet_flux(dataElement);
            end
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Building the right hand side vector
        %------------------------------------------------------------------
        
        function obj = building_Y_vector(obj)
            Nel    = obj.ProblemData.Nel;
            
            index_Elements = obj.ProblemData.inElementVectorLocation;
            index_lines    = obj.ProblemData.lineVectorLocation;
            index_points   = obj.ProblemData.pointVectorLocation;
            elementsData   = obj.ProblemData.elements;
            
            obj.Y.vector = ...
                         zeros(1,obj.ProblemData.pointVectorLocation(end));
            for k = 1:Nel
                EL = obj.Y.CurrentDensity{k}+...
                         obj.Y.Magnetisation{k};
                % Getting the connectivity of the lines and points of
                % k-th element
                if any(k==elementsData.periodic)
                    line_element  = elementsData.periodic_lines(k,:);
                    point_element = elementsData.periodic_points(k,:);
                else
                    line_element  = abs(elementsData.lines(k,:));
                    line_element_no_abs = elementsData.lines(k,:);
                    point_element = elementsData.points(k,:);
                end
               
               % Boundary Conditions
                dir_points = eval(['[',xml2matlab(obj.xmlContent...
                 ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
                dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
                 ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
                % Filling the inElement content
                obj.Y.vector(index_Elements{k}) = int_el(EL,'intElement');
                
                % Filling the line content
                for ii = 1:4
                    if any(line_element(ii)==dir_lines)
                    else
                        tmp_line = ...
                        cell2mat(int_el(EL,ii,'lines'));
                        if line_element_no_abs(ii)<0
                            obj.Y.vector(index_lines{line_element(ii)})=...
                            obj.Y.vector(index_lines{line_element(ii)})+...
                            fliplr(tmp_line);
                            
                        else
                            obj.Y.vector(index_lines{line_element(ii)})=...
                            obj.Y.vector(index_lines{line_element(ii)})+...
                            tmp_line;
                        end
                    end
                end
                
                % Filling the points content
                for ii = 1:4
                    if any(point_element(ii)==dir_points)
                    else
                       obj.Y.vector(index_points(point_element(ii))) = ...
                       obj.Y.vector(index_points(point_element(ii))) + ...
                       int_el(EL,ii,'points');
                    end
                end
            end
            
            % Filling the imposed Flux boundary condition
            
            end
        %==================================================================
        
        %------------------------------------------------------------------
        % Build the global matrix
        %------------------------------------------------------------------
        function obj = global_matrix(obj)
            % Data initialisation
            Nel    = obj.ProblemData.Nel;
            Nl     = obj.ProblemData.Nl;
            Np     = obj.ProblemData.Np;
            
            E_sum  = zeros(obj.ProblemData.pointVectorLocation(end));
            I      = eye(obj.ProblemData.pointVectorLocation(end));
            
            % Boundary Conditions
                dir_points = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
                dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
                per_points = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
                per_lines  = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);
          
            
            for k = 1:Nel
                % data for the element matrix
                dataElement.Lx   = obj.metrics.Lx{k};
                dataElement.Ly   = obj.metrics.Ly{k};
                dataElement.Lxp  = obj.metrics.Lxp{k};
                dataElement.Lyp  = obj.metrics.Lyp{k};
                dataElement.Xcsi = obj.metrics.Xcsi{k};
                dataElement.Ycsi = obj.metrics.Ycsi{k};
                dataElement.Xeta = obj.metrics.Xeta{k};
                dataElement.Yeta = obj.metrics.Yeta{k};
                dataElement.J    = obj.metrics.J{k};
                dataElement.nu   = 1./obj.Materials.Permeability{k};
                dataElement.W    = obj.metrics.W{k};
                dataElement.w_csi = obj.metrics.w_csi{k};
                dataElement.w_eta = obj.metrics.w_eta{k};
                elementsData     = obj.ProblemData.elements;
                
                % Connectivity of elements
                
                if any(k==elementsData.periodic)
                    line_element  = elementsData.periodic_lines(k,:);
                    point_element = elementsData.periodic_points(k,:);
                else
                    line_element  = elementsData.lines(k,:);
                    point_element = elementsData.points(k,:);
                end
                
                % Initalisation of indexes
                tmp        = int_el(dataElement.J,'address');
                yel        = tmp{1};
                yb         = tmp{2};
                yc         = tmp{3};
                i_node     = obj.ProblemData.inElementVectorLocation{k};
                e_node     = obj.ProblemData.lineVectorLocation;
                c_node     = obj.ProblemData.pointVectorLocation;
                %----------------------------------------------------------
                E = element_matrix(dataElement);
                %----------------------------------------------------------
                % Inverse axis
                in_ln = find(line_element<0);
                if in_ln>0
                   inversed_line =  yb{in_ln};
                   LL = E(inversed_line,inversed_line);
                   EL = E(:,inversed_line);
                   LE = E(inversed_line,:);
                   E(inversed_line,:) = flipud(LE);
                   E(:,inversed_line) = fliplr(EL);
                   E(inversed_line,inversed_line) = rot90(LL,2);
                end
         
        %------------------------------------------------------------
        E_sum(i_node,i_node) =  E(yel,yel);
        %=================================================================
        jj = abs(line_element); ll = point_element;
        % inElement Points
        
        E_sum(i_node,c_node(ll)) = E_sum(i_node,c_node(ll)) + E(yel,yc);
        
        % Points inElement
        E_sum(c_node(ll),i_node) = E_sum(c_node(ll),i_node) + E(yc,yel);
        
        % Points Points
        E_sum(c_node(ll),c_node(ll)) = ...
                                   E_sum(c_node(ll),c_node(ll)) + E(yc,yc);
                               
               
                               
            for ii = 1:4
               % element edges
               E_sum(i_node,e_node{jj(ii)}) =   E(yel,yb{ii}); 

               % edges element
               
               E_sum(e_node{jj(ii)},i_node) =   E(yb{ii},yel);
               
               % edges edges
               for n = 1:4
                   E_sum(e_node{jj(ii)},e_node{jj(n)}) = E_sum(...
                       e_node{jj(ii)},e_node{jj(n)}) + E(yb{ii},yb{n});
               end
            
               % edges corners
               E_sum(e_node{jj(ii)},c_node(ll)) = ...
               E_sum(e_node{jj(ii)},c_node(ll)) + E(yb{ii},yc);
        
               % corners edges
               E_sum(c_node(ll),e_node{jj(ii)}) = ...
               E_sum(c_node(ll),e_node{jj(ii)}) + E(yc,yb{ii});
            end

            end

            % dirichlet
            % points
            for k = 1:Np
                if any(k==dir_points)
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==dir_lines)
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end
            
            % periodic
            % points
            for k = 1:Np
                if any(k==per_points)
%                    E_sum(:,c_node(k))  = 0;
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==per_lines)
%                    E_sum(:,e_node{k}) = 0;
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end

            obj.Global_Matrix = E_sum;
        end
        
        %==================================================================
        
        
        %------------------------------------------------------------------
        % Build the global matrix without the permeability
        %------------------------------------------------------------------
        function G = global_matrix_no_nu(obj,nu)
            % Data initialisation
            Nel    = obj.ProblemData.Nel;
            Nl     = obj.ProblemData.Nl;
            Np     = obj.ProblemData.Np;
            
            E_sum  = zeros(obj.ProblemData.pointVectorLocation(end));
            I      = zeros(obj.ProblemData.pointVectorLocation(end));
            
            % Boundary Conditions
                dir_points = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
                dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
                per_points = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
                per_lines  = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);
          
            
            for k = 1:Nel
                % data for the element matrix
                dataElement.Lx   = obj.metrics.Lx{k};
                dataElement.Ly   = obj.metrics.Ly{k};
                dataElement.Lxp  = obj.metrics.Lxp{k};
                dataElement.Lyp  = obj.metrics.Lyp{k};
                dataElement.Xcsi = obj.metrics.Xcsi{k};
                dataElement.Ycsi = obj.metrics.Ycsi{k};
                dataElement.Xeta = obj.metrics.Xeta{k};
                dataElement.Yeta = obj.metrics.Yeta{k};
                dataElement.J    = obj.metrics.J{k};
                dataElement.W    = obj.metrics.W{k};
                dataElement.nu   = nu{k};
                elementsData     = obj.ProblemData.elements;
                
                % Connectivity of elements
                
                if any(k==elementsData.periodic)
                    line_element  = elementsData.periodic_lines(k,:);
                    point_element = elementsData.periodic_points(k,:);
                else
                    line_element  = elementsData.lines(k,:);
                    point_element = elementsData.points(k,:);
                end
                
                % Initalisation of indexes
                tmp        = int_el(dataElement.J,'address');
                yel        = tmp{1};
                yb         = tmp{2};
                yc         = tmp{3};
                i_node     = obj.ProblemData.inElementVectorLocation{k};
                e_node     = obj.ProblemData.lineVectorLocation;
                c_node     = obj.ProblemData.pointVectorLocation;
                %----------------------------------------------------------
                E = element_matrix(dataElement);
                %----------------------------------------------------------
                % Inverse axis
                in_ln = find(line_element<0);
                if in_ln>0
                   inversed_line =  yb{in_ln};
                   LL = E(inversed_line,inversed_line);
                   EL = E(:,inversed_line);
                   LE = E(inversed_line,:);
                   E(inversed_line,:) = flipud(LE);
                   E(:,inversed_line) = fliplr(EL);
                   E(inversed_line,inversed_line) = rot90(LL,2);
                end
         
        %------------------------------------------------------------
        E_sum(i_node,i_node) =  E(yel,yel);
        %=================================================================
        jj = abs(line_element); ll = point_element;
        % inElement Points
        
        E_sum(i_node,c_node(ll)) = E_sum(i_node,c_node(ll)) + E(yel,yc);
        
        % Points inElement
        E_sum(c_node(ll),i_node) = E_sum(c_node(ll),i_node) + E(yc,yel);
        
        % Points Points
        E_sum(c_node(ll),c_node(ll)) = ...
                                   E_sum(c_node(ll),c_node(ll)) + E(yc,yc);
                               
               
                               
            for ii = 1:4
               % element edges
               E_sum(i_node,e_node{jj(ii)}) =   E(yel,yb{ii}); 

               % edges element
               
               E_sum(e_node{jj(ii)},i_node) =   E(yb{ii},yel);
               
               % edges edges
               for n = 1:4
                   E_sum(e_node{jj(ii)},e_node{jj(n)}) = E_sum(...
                       e_node{jj(ii)},e_node{jj(n)}) + E(yb{ii},yb{n});
               end
            
               % edges corners
               E_sum(e_node{jj(ii)},c_node(ll)) = ...
               E_sum(e_node{jj(ii)},c_node(ll)) + E(yb{ii},yc);
        
               % corners edges
               E_sum(c_node(ll),e_node{jj(ii)}) = ...
               E_sum(c_node(ll),e_node{jj(ii)}) + E(yc,yb{ii});
            end

            end

            % dirichlet
            % points
            for k = 1:Np
                if any(k==dir_points)
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==dir_lines)
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end
            
            % periodic
            % points
            for k = 1:Np
                if any(k==per_points)
%                    E_sum(:,c_node(k))  = 0;
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==per_lines)
%                    E_sum(:,e_node{k}) = 0;
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end

            G = E_sum;
        end
        %------------------------------------------------------------------
        % Build the global matrix without the permeability
        %------------------------------------------------------------------
        function G = global_matrix_contour(obj,nu)
            % Data initialisation
            Nel    = obj.ProblemData.Nel;
            Nl     = obj.ProblemData.Nl;
            Np     = obj.ProblemData.Np;
            
            E_sum  = zeros(obj.ProblemData.pointVectorLocation(end));
            I      = zeros(obj.ProblemData.pointVectorLocation(end));
            
            % Boundary Conditions
                dir_points = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
                dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
              ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
                per_points = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
                per_lines  = eval(['[',xml2matlab(obj.xmlContent...
             ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);
          
            
            for k = 1:Nel
                % data for the element matrix
                dataElement.Lx    = obj.metrics.Lx{k};
                dataElement.Ly    = obj.metrics.Ly{k};
                dataElement.Lxp   = obj.metrics.Lxp{k};
                dataElement.Lyp   = obj.metrics.Lyp{k};
                dataElement.Xcsi  = obj.metrics.Xcsi{k};
                dataElement.Ycsi  = obj.metrics.Ycsi{k};
                dataElement.Xeta  = obj.metrics.Xeta{k};
                dataElement.Yeta  = obj.metrics.Yeta{k};
                dataElement.J     = obj.metrics.J{k};
                dataElement.W     = obj.metrics.W{k};
                dataElement.nu    = nu{k};
                elementsData      = obj.ProblemData.elements;
                
                % Connectivity of elements
                
                if any(k==elementsData.periodic)
                    line_element  = elementsData.periodic_lines(k,:);
                    point_element = elementsData.periodic_points(k,:);
                else
                    line_element  = elementsData.lines(k,:);
                    point_element = elementsData.points(k,:);
                end
                
                % Initalisation of indexes
                tmp        = int_el(dataElement.J,'address');
                yel        = tmp{1};
                yb         = tmp{2};
                yc         = tmp{3};
                i_node     = obj.ProblemData.inElementVectorLocation{k};
                e_node     = obj.ProblemData.lineVectorLocation;
                c_node     = obj.ProblemData.pointVectorLocation;
                %----------------------------------------------------------
                E = element_matrix(dataElement);
                %----------------------------------------------------------
                % Inverse axis
                in_ln = find(line_element<0);
                if in_ln>0
                   inversed_line =  yb{in_ln};
                   LL = E(inversed_line,inversed_line);
                   EL = E(:,inversed_line);
                   LE = E(inversed_line,:);
                   E(inversed_line,:) = flipud(LE);
                   E(:,inversed_line) = fliplr(EL);
                   E(inversed_line,inversed_line) = rot90(LL,2);
                end
         
        %------------------------------------------------------------
        E_sum(i_node,i_node) =  E(yel,yel);
        %=================================================================
        jj = abs(line_element); ll = point_element;
        % inElement Points
        
        E_sum(i_node,c_node(ll)) = E_sum(i_node,c_node(ll)) + E(yel,yc);
        
        % Points inElement
        E_sum(c_node(ll),i_node) = E_sum(c_node(ll),i_node) + E(yc,yel);
        
        % Points Points
        E_sum(c_node(ll),c_node(ll)) = ...
                                   E_sum(c_node(ll),c_node(ll)) + E(yc,yc);
                               
               
                               
            for ii = 1:4
               % element edges
               E_sum(i_node,e_node{jj(ii)}) =   E(yel,yb{ii}); 
% 
%                % edges element
%                
               E_sum(e_node{jj(ii)},i_node) =   E(yb{ii},yel);
               
               % edges edges
               for n = 1:4
                   E_sum(e_node{jj(ii)},e_node{jj(n)}) = E_sum(...
                       e_node{jj(ii)},e_node{jj(n)}) + E(yb{ii},yb{n});
               end
            
               % edges corners
               E_sum(e_node{jj(ii)},c_node(ll)) = ...
               E_sum(e_node{jj(ii)},c_node(ll)) + E(yb{ii},yc);
        
               % corners edges
               E_sum(c_node(ll),e_node{jj(ii)}) = ...
               E_sum(c_node(ll),e_node{jj(ii)}) + E(yc,yb{ii});
            end

            end

            % dirichlet
            % points
            for k = 1:Np
                if any(k==dir_points)
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==dir_lines)
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end
            
            % periodic
            % points
            for k = 1:Np
                if any(k==per_points)
%                    E_sum(:,c_node(k))  = 0;
                   E_sum(c_node(k),:)  = I(c_node(k),:);
                end
            end

            % lines
            for k = 1:Nl
                if any(k==per_lines)
%                    E_sum(:,e_node{k}) = 0;
                   E_sum(e_node{k},:) = I(e_node{k},:);
                end
            end

            G = E_sum;
        end
        %------------------------------------------------------------------
        % Update the sources
        %------------------------------------------------------------------
        function obj = updateSources(obj,NewSources)
            obj.Sources     = NewSources;
        end
        
        %------------------------------------------------------------------
        % Update the materials
        %------------------------------------------------------------------
        function obj = updateMaterials(obj,NewMaterials)
            obj.Materials     = NewMaterials;
        end

    end
    
end


function out = xml2matlab(varargin)
         mode          = varargin{nargin};
         xmlElement    = varargin{1};
         ElementName   = varargin{2};
         ElementNumber = varargin{3};
switch mode
    % Extracts all the attributes for certain element and returns the name
    % and its value in a 2 column cell
    case 'Attributes'
        Attributes = xmlElement.getElementsByTagName(ElementName). ...
            item(ElementNumber).getAttributes;
        AttributeNumber = Attributes.getLength;
        Parameters = cell(AttributeNumber,2);
        for k = 1:AttributeNumber
            Parameters{k,1} = char(Attributes.item(k-1).getName);
            Parameters{k,2} = char(Attributes.item(k-1).getValue);
        end
        out = Parameters;
    % Extracts only selected attribute and returns its value
    case 'Attribute'
        AttributeName = varargin{4};
        Attribute = char(xmlElement.getElementsByTagName(ElementName). ...
            item(ElementNumber).getAttribute(AttributeName));
        out = Attribute;
        
end
end

%------------------------------------------------------------------
% Element Matrix implementation
%------------------------------------------------------------------
function Eij = element_matrix(dataElement)
Lx   = dataElement.Lx;    Ly   = dataElement.Ly;
Lxp  = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta = dataElement.Xeta;  Yeta = dataElement.Yeta;
J    = dataElement.J;
nu   = dataElement.nu;
W = dataElement.W;


A = nu.*(Yeta.^2 + Xeta.^2)./J;         A  = diag(A(:));
B = nu.*(Yeta.*Ycsi + Xeta.*Xcsi)./J;   B  = diag(B(:));
C = nu.*(Ycsi.^2 + Xcsi.^2)./J;         C  = diag(C(:));

W = diag(W(:));


Eij =   - (Lxp*W*(A*Lx - B*Ly) + Lyp*W*(C*Ly - B*Lx));

end 


function Eij = magnet_flux(dataElement)

Lxp   = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi  = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta  = dataElement.Xeta;  Yeta = dataElement.Yeta;
nu    = dataElement.nu;
w_csi = dataElement.w_csi;
w_eta = dataElement.w_eta;

F1c   =   Yeta.*(-dataElement.My).*nu - Xeta.*dataElement.Mx.*nu;
F2c   = - Ycsi.*(-dataElement.My).*nu + Xcsi.*dataElement.Mx.*nu;

W    = w_eta*w_csi;

Eij =  -(Lxp*diag(W(:))*(F1c(:)) + Lyp*diag(W(:))*F2c(:));
Eij = reshape(Eij,size(nu));
end