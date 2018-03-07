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
               obj.Fourier = FourierElement(obj.Physics.xmlContent,obj.Geometry);
                % Building the Fourier to space tranformation matrix
                Espace = obj.Fourier.fourier_space_matrix...
                                                 (obj.Problem.ProblemData);
                
                Efrequency = obj.Fourier.fourier_frequency_matrix...
                                                 (obj.Problem.ProblemData);
                % Storing matrices
                obj.Fourier_matrix.Espace        = Espace;
                obj.Fourier_matrix.Efrequency    = Efrequency;
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
        obj = interpolate_on_FEM(obj,filenames,no_elements)
    end
    
end
