%==========================================================================
% Problem.m
% Created: 27.10.2017 - 14:45:01
% By: M. Curti
%
% Building the global matrix and the right hand side source vector
%==========================================================================
classdef Problem < matlab.mixin.SetGet
    %PROBLEM Builds the matrices for the problem
    %   Detailed explanation goes here
    
    properties
        Y;             % Right hand side of the equation
        ProblemData;   % Various parameters of the problem
        Sources;       % Sources in the geometry
        metrics        % Metrics of the geometry
        mappings        % Metrics of the geometry
        xmlContent     % Info from the xml file
        Materials      % Materials of the problem
        tmp_Element_Matrix
        Global_Matrix
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = Problem(G, Ph)
            obj.ProblemData = Ph.ProblemData;
            obj.Sources     = Ph.Sources;
            obj.metrics     = G.metrics;
            obj.mappings    = G.mappings;
            obj.xmlContent  = G.GeometryElement;
            obj.Materials   = Ph.Materials;
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
                    line_element  = abs(elementsData.periodic_lines(k,:));
                    point_element = elementsData.periodic_points(k,:);
                    line_element_no_abs = line_element;
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
        obj = global_matrix(obj)
        
        %------------------------------------------------------------------
        % Build the global matrix in axysimmetric coordinates
        %------------------------------------------------------------------
        obj = global_matrix_axi(obj)
        %------------------------------------------------------------------
        % Build the global matrix without the permeability
        %------------------------------------------------------------------
        JacMatrix = global_matrix_jacobian(obj,kk,Potential)
        %------------------------------------------------------------------
        % Build the global matrix without the permeability
        %------------------------------------------------------------------
        [G, J] = global_matrix_and_jacobian(obj,kk,Potential)
        %------------------------------------------------------------------
        % Build the global matrix without the permeability
        %------------------------------------------------------------------
        obj = global_matrix_contour(obj,nu)
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
        
        %------------------------------------------------------------------
        % jacobian and quadrature product
        %------------------------------------------------------------------
        qj = quad_jacobian(obj)
        
        %------------------------------------------------------------------
        % inverse r
        %------------------------------------------------------------------
        rm1 = inverse_r(obj)

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




function Eij = magnet_flux(dataElement)

Lxp   = dataElement.Lxp;   Lyp  = dataElement.Lyp;
Xcsi  = dataElement.Xcsi;  Ycsi = dataElement.Ycsi;
Xeta  = dataElement.Xeta;  Yeta = dataElement.Yeta;
nu    = dataElement.nu;
w_csi = dataElement.w_csi;
w_eta = dataElement.w_eta;

F1c   =   Yeta.*(-dataElement.My).*nu - Xeta.*dataElement.Mx.*nu;
F2c   = - Ycsi.*(-dataElement.My).*nu + Xcsi.*dataElement.Mx.*nu;

W     = w_eta*w_csi;

Eij =  -(Lxp*diag(W(:))*(F1c(:)) + Lyp*diag(W(:))*F2c(:));
Eij = reshape(Eij,size(nu));
end