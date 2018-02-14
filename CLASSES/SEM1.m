%==========================================================================
% SEM.m
% Created: 29.10.2017 - 22:05:25
% By: M. Curti
%
% The main class of the method
%==========================================================================
classdef SEM < matlab.mixin.SetGet
    %SEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Geometry; % The instance of the Geometry class
        Physics;  % The instance of the Physics class
        Problem;  % The instance of the Problem class
        PProcessing; % The instance of the post processing class
        Fourier;     % The Fourier class
        Fourier_matrix; % The set of fourier matrices
        NonlinearSolver; % The settings for the nonlinear solver
        SEM_interpolation;
    end
    
    methods
        
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = SEM(filename)
            % The constructor starts with reading the XML file where all
            % SEM problem info is available
            
            % Loading the default settings
            % - Settings for the nonlinear solver
            obj.NonlinearSolver.max_iter = 30;   % Maximum iterations for
                                               % the nonlinear solver
            obj.NonlinearSolver.tol    = 1e-4;   % Maximum tolerance for
                                               % the nonlinear solver
            obj.NonlinearSolver.display = 'off'; % Show or hide the solver 
                                               % iterations
            % Storing the instance of the Geometry class
            disp('SEM - Building the geometry')
            obj.Geometry = Geometry(filename);
            
            % Storing the instance of the Physics class
            disp('SEM - Loading the physics data')
            obj.Physics = Physics(obj.Geometry);
            
            obj.Physics = obj.Physics.load_magnetic_materials;
            obj.Physics = obj.Physics.load_current_density;
            
            % Storing the instance of the Problem class
            disp('SEM - Preparing the data for the linear problem')
            dataProblem.ProblemData = obj.Physics.ProblemData;
            dataProblem.Sources     = obj.Physics.Sources;
            dataProblem.metrics     = obj.Geometry.metrics;
            dataProblem.xmlContent  = obj.Geometry.GeometryElement;
            dataProblem.Materials   = obj.Physics.Materials;
            
            obj.Problem = Problem(dataProblem);
            
            % Load the sources and build the Y vector
            disp('SEM - Load the sources and build the rhs vector')
            obj.Problem = obj.Problem.load_Y_sources;
            obj.Problem = obj.Problem.building_Y_vector;
            
            %==============================================================
            % Fourier Regions
            %==============================================================
            try
                % Geting general parameters for Fourier
                Q     = eval(xml2matlab(obj.Problem.xmlContent...
                    ,'Fourier',0,'Harmonics','Attribute'));
                type  = xml2matlab(obj.Problem.xmlContent...
                    ,'Fourier',0,'type','Attribute');
                start = eval(xml2matlab(obj.Problem.xmlContent...
                    ,'Fourier',0,'start','Attribute'));
                tau   = eval(xml2matlab(obj.Problem.xmlContent...
                    ,'Fourier',0,'tau','Attribute'));
                fel   = eval(xml2matlab(obj.Problem.xmlContent...
                    ,'Fourier',0,'fel','Attribute'));
                
                % Geting parameters per domain
                parameter_loader(obj.Geometry.parameters)
                domain.heights    = zeros(fel,2);
                domain.points_top = cell(fel);
                domain.lines_top  = cell(fel);
                
                for k = 1:fel
                    domain.heights(k,:) = ...
                        eval(['[' xml2matlab(obj.Problem.xmlContent...
                        ,'domain',k-1,'heights','Attribute') '];']);
                    domain.lines_top{k} = ...
                        eval(['[' xml2matlab(obj.Problem.xmlContent...
                        ,'domain',k-1,'lines_top','Attribute') '];']);
                    domain.elements_top{k} = ...
                        eval(['[' xml2matlab(obj.Problem.xmlContent...
                        ,'domain',k-1,'elements_top','Attribute') '];']);
                end
                
                % Formating the data for Fourier class
                Elements           = struct;
                Elements.x_start   = start;
                Elements.tau       = tau;
                Elements.Harmonics = Q;
                Elements.heights   = domain.heights;
                Elements.type      = type;
                obj.Fourier        = fourierElements(Elements);
                Geom_data.lines    = obj.Geometry.lines;
                Geom_data.xi       = obj.Geometry.xi;
                Geom_data.metrics  = obj.Geometry.metrics;
                %
                FourierData.connectedLines = domain.lines_top{1};
                FourierData.connectedElements = domain.elements_top{1};
                % Building the Fourier to space tranformation matrix
                [Espace, in_line] = obj.Fourier.fourier_space_matrix...
                           (obj.Problem.ProblemData,Geom_data,FourierData);
                
                Efrequency = obj.Fourier.fourier_frequency_matrix...
                           (obj.Problem.ProblemData,Geom_data,FourierData);
                
                       
                % Storing matrices
                obj.Fourier_matrix.Espace     = Espace;
                obj.Fourier_matrix.in_line    = in_line;
                obj.Fourier_matrix.Efrequency = Efrequency;
            catch
                disp('No Fourier domains detected')
            end
        end
        
        %------------------------------------------------------------------
        % Solve function; Bulding the right and the left hand side of the 
        % linear problem and solve the problem
        %------------------------------------------------------------------
        obj = solve(obj)
        %------------------------------------------------------------------
        % Function to interpolate on FEM nodes, used to check the
        % correctness of the problem
        %------------------------------------------------------------------
        function obj = interpolate_on_FEM(obj,filenames,no_elements)
            
            FluxData  = Flux2DmeshData(filenames);
            FEM       = FluxData.FEMfilesInOne;
            try
                FEM.noelements = no_elements;
            catch
            end
            sem = 0;
            fileName = sprintf('lxy%d.mat',numel(FEM.all_x_nodes));
            if exist(fileName,'file')
                load(fileName)
                FEM.int  = sem.int;
                obj.SEM_interpolation.solutions = sem;
                obj.SEM_interpolation.FEM       = FEM.all_values;
                obj.SEM_interpolation.TRI = FEM.TRI;
                obj.SEM_interpolation.all_x_nodes = FEM.all_x_nodes;
                obj.SEM_interpolation.all_y_nodes = FEM.all_y_nodes;
            else
                sem = obj.PProcessing.InterpolateOnFEMNodes(FEM);
                
                obj.SEM_interpolation.solutions = sem;
                obj.SEM_interpolation.FEM       = FEM.all_values;
                obj.SEM_interpolation.TRI = FEM.TRI;
                obj.SEM_interpolation.all_x_nodes = FEM.all_x_nodes;
                obj.SEM_interpolation.all_y_nodes = FEM.all_y_nodes;
                save(fileName,'sem');
            end
            
        end
    end
    
end
