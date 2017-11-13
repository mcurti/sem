%==========================================================================
% SEM.m
% Created: 29.10.2017 - 22:05:25
% By: M. Curti
%
% The main class of the method
%==========================================================================
classdef SEM
    %SEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Geometry; % The instance of the Geometry class
        Physics;  % The instance of the Physics class
        Problem;  % The instance of the Problem class
        PProcessing; % The instance of the post processing class
        NonlinearSolver; % The settings for the nonlinear solver
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
            obj.NonlinearSolver.max_iter = 30; % Maximum iterations for the
                                               % nonlinear solver
            obj.NonlinearSolver.tol    = 1e-6; % Maximum tolerance for the 
                                               % nonlinear solver
            
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
        end
        
        %------------------------------------------------------------------
        % Solve function; Bulding the right and the left hand side of the 
        % linear problem and solve the problem
        %------------------------------------------------------------------
        function obj = solve(obj)
            
            % Initiating the variables
            
            iter = obj.NonlinearSolver.max_iter;
            
            err = zeros(1,iter);
            
            Y_mag = ...
                zeros(obj.Problem.ProblemData.pointVectorLocation(end),1);
            prev_PHI = Y_mag;
            
            
            dataPostProcessing.xmlContent   = ...
                obj.Geometry.GeometryElement;
            dataPostProcessing.ProblemData  = obj.Physics.ProblemData;
            dataPostProcessing.mappings     = obj.Geometry.mappings;
            dataPostProcessing.metrics      = obj.Geometry.metrics;
            dataPostProcessing.lines        = obj.Geometry.lines;
            dataPostProcessing.points       = obj.Geometry.points;
            dataPostProcessing.xi           = obj.Geometry.xi;
            dataPostProcessing.elements     = obj.Geometry.elements;
            dataPostProcessing.Permeability = ...
                obj.Problem.Materials;
            
            % Load the sources and build the Y vector
            disp('SEM - Load the sources and build the rhs vector')
            obj.Problem = obj.Problem.load_Y_sources;
            obj.Problem = obj.Problem.building_Y_vector;
            
            % Newton Raphson solver
            tic
            fprintf('Starting the nonlinear solver %.4f \n',toc);
            
            for ii = 1:iter
                
                fprintf('Building the global matrix at iteration %d \n'...
                    ,ii);
                obj.Problem = obj.Problem.global_matrix;
                
                fprintf('Start solving the linear system at time %.4f \n'...
                    ,toc);
                s_time = toc;
                PHI = obj.Problem.Global_Matrix\...
                                             (obj.Problem.Y.vector'+Y_mag);
                fprintf('Linear system solved in  %.4f seconds \n', ...
                    toc - s_time);
                %% pause
                
                
                dataPostProcessing.PHI          = PHI;
                
                disp('Evaluating the convergece')
                
                err(ii) = max(abs(PHI-prev_PHI));
                prev_PHI = PHI;
                
                obj.PProcessing = PostProcessing(dataPostProcessing);
                obj.PProcessing = obj.PProcessing.compute_B;
                
                [MU, KK]  = obj.PProcessing.get_next_permeability(...
                    obj.Problem.Materials,obj.Geometry.GeometryElement);
                
                
                new_material.Permeability = MU;
                obj.Problem = obj.Problem.updateMaterials(new_material);
                MagMatrix = obj.Problem.global_matrix_contour(KK);
                Y_mag = MagMatrix*PHI;
                
                if err(ii) < obj.NonlinearSolver.tol
                    disp('The solution converged to desired tolerance')
                    break
                end
                
            end
            
%             figure(1)
%             clf
%             semilogy(1:ii,err(1:ii))
%             pause
            
            if max(abs(Y_mag)>0)
               PHI_rem = obj.Problem.Global_Matrix\Y_mag;
            else
                PHI_rem = PHI*0;
            end
                dataPostProcessing.PHI = PHI_rem;
                
            Pp = PostProcessing(dataPostProcessing);
            
            set(obj.PProcessing,'remPotential',Pp.Potential);
                
        end
    end
    
end

