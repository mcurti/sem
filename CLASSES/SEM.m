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
        end
        
        %------------------------------------------------------------------
        % Solve function; Bulding the right and the left hand side of the 
        % linear problem and solve the problem
        %------------------------------------------------------------------
        function obj = solve(obj)
            
            % Initiating the variables
            
            iter = obj.NonlinearSolver.max_iter;
            Nel  = obj.Problem.ProblemData.Nel;
            err = zeros(1,iter);
            obj.NonlinearSolver.MU = cell(iter,Nel);
            obj.NonlinearSolver.K  = cell(iter,Nel);
            
            Y_mag = ...
                zeros(obj.Problem.ProblemData.pointVectorLocation(end),1);
%             prev_PHI = Y_mag;
            
            
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
            
            err_el = zeros(1,Nel);
            init_phi = cell(1,Nel);
            % Count the periodic unknowns
            per_points = eval(['[',xml2matlab(obj.Problem.xmlContent...
             ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
                per_lines  = eval(['[',xml2matlab(obj.Problem.xmlContent...
             ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);
         
            matrix_index = true(...
                obj.Physics.ProblemData.inElementsUnknowns + ...
                obj.Physics.ProblemData.linesUnknowns      + ...
                obj.Physics.ProblemData.Np,1);
            e_node = obj.Physics.ProblemData.lineVectorLocation;
            c_node = obj.Physics.ProblemData.pointVectorLocation;
            
            for k = per_points
                
                matrix_index(c_node(k),:)  = false;
            end

            % lines
            for k = per_lines
                
                matrix_index(e_node{k},:) = false;
            end
            
            for k = 1:Nel
                init_phi{k} = zeros(size(obj.Problem.metrics.J{k}));
                
            end
            PHI = zeros(size(matrix_index));
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
                S = sparse(obj.Problem.Global_Matrix);
                
                PHI(matrix_index) = S(matrix_index,matrix_index)...
                    \(obj.Problem.Y.vector(matrix_index)' + ...
                     Y_mag(matrix_index));
                 
                fprintf('Linear system solved in  %.4f seconds \n', ...
                    toc - s_time);
                %% pause
                
                dataPostProcessing.PHI          = PHI;
                
                disp('Evaluating the convergece')
                
                
                obj.PProcessing = PostProcessing(dataPostProcessing);
                obj.PProcessing = obj.PProcessing.compute_B;
                
                [MU, KK, ne]  = obj.PProcessing.get_next_permeability(...
                    obj.Problem.Materials,obj.Geometry.GeometryElement);
                
                
                for  k = 1:Nel
                    err_el(k) = sqrt(sum(((init_phi{k}(:) - obj.PProcessing.Potential{k}(:)).^2).*obj.Geometry.metrics.J{k}(:).*obj.Geometry.metrics.W{k}(:)));
                    obj.NonlinearSolver.MU{ii,k} = MU{k};
                    obj.NonlinearSolver.K{ii,k}  = KK{k};
                end
                
                if any(ii==[5 10 15])
                    for k = 1:Nel
                        MU{k} = MU{k} - (-obj.NonlinearSolver.MU{ii-1,k} + MU{k})*0.1;
                        KK{k} = KK{k} - (-obj.NonlinearSolver.K{ii-1,k}  + KK{k})*0.1;
                    end
                end
                err(ii) = max(err_el);
                
%                 prev_PHI = PHI;
                init_phi = obj.PProcessing.Potential;
                %----------------------------------------------------------
                % Display the iterations
                %----------------------------------------------------------
                if strcmp(obj.NonlinearSolver.display,'off')
                elseif strcmp(obj.NonlinearSolver.display,'on')
                    figure(100)
                    clf
                    semilogy(1:ii,err(1:ii))
                    hold on
                    semilogy(1:iter,...
                        obj.NonlinearSolver.tol*ones(1,iter),'-r')
                    hold off
                    xlim([1 iter])
                    ylim([obj.NonlinearSolver.tol*.5 err(1)])
                    xlabel('iterations')
                    ylabel('norm')
                    figure_config(100,15,5,8)
                    drawnow 
                end
                
                %----------------------------------------------------------
                % Check if the system is linear
                %----------------------------------------------------------
                if sum(ne)==0
                    disp('linear system detected, solved with one iteration')
                    break
                else
                    %------------------------------------------------------
                    % Update the permeability
                    %------------------------------------------------------
                    new_material.Permeability = MU;
                    obj.Problem = obj.Problem.updateMaterials(new_material);
                    MagMatrix = obj.Problem.global_matrix_contour(KK);
                    Y_mag = MagMatrix*PHI;
                    
                    %------------------------------------------------------
                    % Check if converged
                    %------------------------------------------------------
                    if err(ii) < obj.NonlinearSolver.tol
                        disp('The solution converged to desired tolerance')
                        break
                    end
                end
                
                %----------------------------------------------------------
                % Display "Not converged" message
                %----------------------------------------------------------
                if iter == ii && err(ii)> obj.NonlinearSolver.tol
                    disp('max iteration achieved')
                    fprintf('tol = %.4g, desired tol = %.4g\n',...
                        err(ii),obj.NonlinearSolver.tol)
                end
                
            end
                        
            if max(abs(Y_mag)>0)
               PHI_rem = obj.Problem.Global_Matrix\Y_mag;
            else
                PHI_rem = PHI*0;
            end
                dataPostProcessing.PHI = PHI_rem;
                
            Pp = PostProcessing(dataPostProcessing);
            
            set(obj.PProcessing,'remPotential',Pp.Potential);
            obj.NonlinearSolver.err = err;
        end
        
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