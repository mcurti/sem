%==========================================================================
% solve.m
% Created: 09.02.2018 - 18:26:58
% By: M. Curti
%
% the set of linearly independent functions for polar coordinate system
%==========================================================================
function obj = solve(obj)

% Initiating the variables

iter = obj.NonlinearSolver.max_iter;
Nel  = obj.Problem.ProblemData.Nel;
err = zeros(1,iter);
obj.NonlinearSolver.MU = cell(iter,Nel);
obj.NonlinearSolver.K  = cell(iter,Nel);

if isempty(obj.Fourier)
    Y_mag = ...
        zeros(obj.Problem.ProblemData.pointVectorLocation(end),1);
else
    Y_mag = [zeros(...
        obj.Problem.ProblemData.pointVectorLocation(end),1); ...
        zeros(size(obj.Fourier_matrix.Efrequency,1),1)];
end
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


tot_unknowns = (obj.Physics.ProblemData.inElementsUnknowns + ...
    obj.Physics.ProblemData.linesUnknowns + obj.Physics.ProblemData.Np);
e_node = obj.Physics.ProblemData.lineVectorLocation;
c_node = obj.Physics.ProblemData.pointVectorLocation;

SEM_index = 1:tot_unknowns;

for k = per_points
    SEM_index(c_node(k),:)  = 0;
end

% lines
for k = per_lines
    
    SEM_index(e_node{k},:) = 0;
end

SEM_index(SEM_index==0) = [];


for k = 1:Nel
    init_phi{k} = zeros(size(obj.Problem.metrics.J{k}));
    
end
% Newton Raphson solver
tic
fprintf('Starting the nonlinear solver %.4f \n',toc);

for ii = 1:iter
    
    fprintf('Building the global matrix at iteration %d \n'...
        ,ii);
    obj.Problem = obj.Problem.global_matrix;
    E = obj.Problem.Global_Matrix;
    Y = obj.Problem.Y.vector;
    % Loading the Fourier matrices
    if isempty(obj.Fourier)
        fourier_index = [];
    else
        space2freq = obj.Fourier_matrix.Efrequency;
        freq2space = obj.Fourier_matrix.Espace;
        fourier_index = (1:(size(freq2space,2))) + size(freq2space,1);
        E = [obj.Problem.Global_Matrix freq2space; space2freq];
        Y = [obj.Problem.Y.vector zeros(1,size(space2freq,1))];
        
    end
    index = [SEM_index fourier_index ];
    PHI = zeros(1,(tot_unknowns+numel(fourier_index)));
    
    fprintf('Start solving the linear system at time %.4f \n'...
        ,toc);
    s_time = toc;
    S = sparse(E);
    
    
    PHI(index) = S(index,index)\(Y(index)' + Y_mag(index));
    
    fprintf('Linear system solved in  %.4f seconds \n', ...
        toc - s_time);
    %% pause
    
    dataPostProcessing.PHI          = PHI;
    dataPostProcessing.PHI(fourier_index) = [];
    
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

%----------------------------------------------------------
% Update the Fourier data
%----------------------------------------------------------
if isempty(obj.Fourier)
else
    Q    = obj.Fourier.Edata.Harmonics;
    Nfel = obj.Fourier.Edata.Nfel;
    X = PHI(fourier_index);
    
    c1 = zeros(1,Q*Nfel); c2 = zeros(1,Q*Nfel); c3 = zeros(1,Q*Nfel);
    c4 = zeros(1,Q*Nfel); %Bx0 = zeros(1,1);
    for k = 1:Nfel
        c1((1:Q) + (k-1)*Q) = X((1:Q) + 4*Q*(k-1));
        c2((1:Q) + (k-1)*Q) = X((1:Q) + Q*(1 + 4*(k-1)));
        c3((1:Q) + (k-1)*Q) = X((1:Q) + Q*(2 + 4*(k-1)));
        c4((1:Q) + (k-1)*Q) = X((1:Q) + Q*(3 + 4*(k-1)));
        %     Bx0(k) = X(end - (k-1));
    end
    
    Bx0 = X(end-1);
    
    Az0 = X(end);
    
    s  = zeros(1,numel(fourier_index));
    obj.Fourier = obj.Fourier.update_coefficients(s, c1, c2, c3, c4, Az0, Bx0);
end
end

