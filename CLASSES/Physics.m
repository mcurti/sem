%==========================================================================
% Physics.m
% Created: 27.10.2017 - 14:42:42
% By: M. Curti
%
% Physics class is loading the physics information from the xml file and
% stores it in the cells, so that the Problem class will get processed
% information. Moreover, it also has methods to update the received
% information.
%==========================================================================
classdef Physics
    %PHYSICS Bulding of the problem for the given geometry.
    %   As the input, it uses the class GEOMETRY
    
    properties
        parameters;    % physical parameters
        G_class;       % Geometry Class
        xmlContent;    % instances from the xml file
        Sources;       % Sources in the geometry
        Materials;     % Materials in the geometry
        ProblemData    % Various data of the problem
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = Physics(G_Class)
            obj.xmlContent = G_Class.GeometryElement;
            obj.G_class    = G_Class;
            % Obtaining the parameters
            obj.parameters = xml2matlab(obj.xmlContent,...
                                     'PhysicalParameters',0,'Attributes');
                                 
            
            % Calculating the problem dimension
            Nel = element_number(obj.xmlContent);
            Nl  = line_number(obj.xmlContent);
            Np  = point_number(obj.xmlContent);
            
            obj.ProblemData.Nel = Nel;
            obj.ProblemData.Nl  = Nl;
            obj.ProblemData.Np  = Np;
            
            obj.ProblemData.ElementSize             = cell(1,Nel);
            obj.ProblemData.inElementVectorLocation = cell(1,Nel);
            obj.ProblemData.lineVectorLocation      = cell(1,Nl);
            obj.ProblemData.GridPoints              = 0;
            obj.ProblemData.elements                = G_Class.elements;
            obj.ProblemData.physical_parameters     = obj.parameters;
            u_tmp = 0;
            for k = 1:Nel
                [M, N] = size(obj.G_class.mappings.Xm{k});
                obj.ProblemData.ElementSize{k} = [M N];
                obj.ProblemData.GridPoints = obj.ProblemData.GridPoints...
                                                  + M*N;
                if (N-2)*(M-2) == 0
                   obj.ProblemData.inElementVectorLocation{k} = [];
                else
                   obj.ProblemData.inElementVectorLocation{k} = ...
                                                   (1:(N-2)*(M-2)) + u_tmp;
                end
                u_tmp = u_tmp + (N-2)*(M-2);
            end
            obj.ProblemData.inElementsUnknowns = u_tmp;
            obj.ProblemData.linesUnknowns = sum(obj.G_class.xi.line_N-1);
            obj.ProblemData.pointsUnknowns = point_number(obj.xmlContent);
            
            for k = 1:Nl
                obj.ProblemData.lineVectorLocation{k} = ...
                    (1:obj.G_class.xi.line_N(k)-1) + u_tmp;
                u_tmp = u_tmp + obj.G_class.xi.line_N(k)-1;
            end
            obj.ProblemData.pointVectorLocation = (1:Np) + u_tmp;
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Load Magnetic Materials
        %------------------------------------------------------------------
        function obj = load_magnetic_materials(obj)
            % Load the number of regions and elements
            Nr = eval(xml2matlab(obj.xmlContent,'Regions'...
                                                     ,0,'nR','Attribute'));
            Nel = element_number(obj.xmlContent);
            Nparam = length(obj.parameters);
            % Loading in the local workspace
            
            check_list = false(1,Nparam);
            
            for k = 1:Nparam
                try
                    eval([obj.parameters{k,1}, '=', ...
                                               obj.parameters{k,2}, ';']);
                catch
                    check_list(k) = true;
                end
            end
            check_list = find(check_list==true);
            % running the parameters with errors
            
            for k = check_list
                    eval([obj.parameters{k,1}, '=', ...
                                               obj.parameters{k,2}, ';']);
            end
            % Updating the default magnetic materials
            obj.Materials.Permeability = cell(1,Nel);
            for k = 1:Nel
                [M, N] = size(obj.G_class.mappings.Xm{k});
                obj.Materials.Permeability{k} = ones(M,N);
            end
            % Updating the assigned properties
            for k = 1:Nr
                ElementList  = eval(['[' xml2matlab(obj.xmlContent,...
                    'region',k-1,'ElementList','Attribute') ']']);
                try
                    permeability = eval(xml2matlab(obj.xmlContent,...
                        'region',k-1,'MagneticMaterial','Attribute'));
                catch
                    
                    permeability = eval(xml2matlab(obj.xmlContent,...
                        'region',k-1,'InitialPermeability','Attribute'));
                end
                for ii = 1:length(ElementList)
                    obj.Materials.Permeability{ElementList(ii)} = ...
                    ones(size(obj.G_class.mappings.Xm{ElementList(ii)}))...
                    *permeability;
                end
            end
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Load Current Density
        %------------------------------------------------------------------
        function obj = load_current_density(obj)
            % Load the number of regions and elements
            Nr = eval(xml2matlab(obj.xmlContent,'Regions'...
                                                     ,0,'nR','Attribute'));
            Nel = element_number(obj.xmlContent);
            % Loading parameters in the local workspace
            Nparam = length(obj.parameters);
            
            check_list = false(1,Nparam);
            % Initialisation of the parameters
            obj.Sources.Magnetisation = cell(2,Nel);
            
            
            for k = 1:Nparam
                try
                    eval([obj.parameters{k,1}, '=', ...
                                               obj.parameters{k,2}, ';']);
                catch
                    check_list(k) = true;
                end
            end
            check_list = find(check_list==true);
            % running the parameters with errors
            
            for k = check_list
                    eval([obj.parameters{k,1}, '=', ...
                                               obj.parameters{k,2}, ';']);
            end
            
            % Updating the default current density
            obj.Sources.CurrentDensity = cell(1,Nel);
            for k = 1:Nel
                obj.Sources.CurrentDensity{k} = ...
                    zeros(size(obj.G_class.mappings.Xm{k}));
                obj.Sources.Magnetisation{1,k} = ...
                    zeros(size(obj.G_class.mappings.Xm{k}));
                obj.Sources.Magnetisation{2,k} = ...
                    zeros(size(obj.G_class.mappings.Xm{k}));
            end
            
            % Updating the properties as assigned
            for k = 1:Nr
                ElementList  = eval(['[' xml2matlab(obj.xmlContent,...
                    'region',k-1,'ElementList','Attribute') ']' ]);
                current = eval(xml2matlab(obj.xmlContent,'region'...
                                     ,k-1,'CurrentDensity','Attribute'));
                magnet  = eval(['[' xml2matlab(obj.xmlContent,'region'...
                                 ,k-1,'Magnetisation','Attribute') ']' ]);
                for ii = 1:length(ElementList)
                    obj.Sources.CurrentDensity{ElementList(ii)} = ...
                    ones(size(obj.G_class.mappings.Xm{ElementList(ii)}))...
                    *current;
                    if isempty(magnet)==false
                        % x - component of the magnetisation
                        
                        obj.Sources.Magnetisation{1,ElementList(ii)} = ...
                            ones(size(obj.G_class.mappings.Xm{...
                            ElementList(ii)}))*magnet(1);
                        % y - component of the magnetisation
                        
                        obj.Sources.Magnetisation{2,ElementList(ii)} = ...
                            ones(size(obj.G_class.mappings.Xm{...
                            ElementList(ii)}))*magnet(2);                        
                    end
                end
            end
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Update Parameters
        %------------------------------------------------------------------
        function obj = setParameter(obj,Parameter,Value)
            obj = change_parameter( obj,Parameter,Value );
        end
    end
    
end


% Loading the number of elements
function Nl = element_number(xml)
    Nl = eval(xml2matlab(xml,'Elements'...
                                                     ,0,'el','Attribute'));
end

function Nl = line_number(xml)
    Nl = eval(xml2matlab(xml,'Lines'...
                                                     ,0,'l','Attribute'));
end

function Nl = point_number(xml)
    Nl = eval(xml2matlab(xml,'Points'...
                                                     ,0,'p','Attribute'));
end