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

tot_unknowns = (obj.Physics.ProblemData.inElementsUnknowns + ...
    obj.Physics.ProblemData.linesUnknowns + obj.Physics.ProblemData.Np);
if isempty(obj.Fourier)
    fourier_index = [];
    spc = tot_unknowns;
    Y_mag = zeros(obj.Problem.ProblemData.pointVectorLocation(end),1);
else
    space2freq = obj.Fourier_matrix.Efrequency;
    freq2space = obj.Fourier_matrix.Espace;
    fourier_index = (1:(size(freq2space,2))) + size(freq2space,1);
    [spr, spc] = size(space2freq); [frr, frc] = size(freq2space);
    Y_mag = [zeros(...
        obj.Problem.ProblemData.pointVectorLocation(end),1); ...
        zeros(size(obj.Fourier_matrix.Efrequency,1),1)];
end
SEM_index = 1:tot_unknowns;

prev_PHI = Y_mag;
% delta_PHI = prev_PHI(SEM_index);

E      = zeros(spc);
JacFin = E;
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

% err_el = zeros(1,Nel);
init_phi = cell(1,Nel);
% Count the periodic unknowns
% per_points = eval(['[',xml2matlab(obj.Problem.xmlContent...
%     ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
% per_lines  = eval(['[',xml2matlab(obj.Problem.xmlContent...
%     ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);


% e_node = obj.Physics.ProblemData.lineVectorLocation;
% c_node = obj.Physics.ProblemData.pointVectorLocation;


% for k = per_points
%     SEM_index(c_node(k),:)  = 0;
% end
%
% % lines
% for k = per_lines
%
%     SEM_index(e_node{k},:) = 0;
% end
%
% SEM_index(SEM_index==0) = [];
obj.Problem = obj.Problem.load_Y_sources;
obj.Problem = obj.Problem.building_Y_vector;

for k = 1:Nel
    init_phi{k} = zeros(size(obj.Problem.metrics.J{k}));

end
% Newton Raphson solver
tic
fprintf('Starting the nonlinear solver %.4f \n',toc);

for ii = 1:iter
    
    % Loading the Fourier matrices
    if ii==1
        %-----------------------------------------------------
        filename = ['PHI' num2str(numel(prev_PHI)) '.mat'];
        try
            load(filename,'PHI')
            
        catch
            disp('No initial solution saved')
            PHI = prev_PHI;
        end
        %-----------------------------------------------------
    end
    %         E = zeros(fourier_index(end));
    
    fprintf('Building the global matrix at iteration %d \n'...
        ,ii);
    dataPostProcessing.PHI          =  PHI(SEM_index);% + 0*delta_PHI;
    %         dataPostProcessing.PHI = [];
    obj.PProcessing = PostProcessing(dataPostProcessing);
    obj.PProcessing = obj.PProcessing.compute_B;
    
%     bh = bh_class(2);
%     figure(11)
%     if ii ==1
%     clf
%     plot(0:.01:2.2,ppval(bh.Extrap_Spline_B,0:.01:2.2)-5.500778160865974e+04)
%     end
%     hold on
%     plot(obj.PProcessing.Flux.abs{1,17}(12,3)*1e3,ppval(bh.Extrap_Spline_B,obj.PProcessing.Flux.abs{1,17}(12,3)*1e3)-5.500778160865974e+04,'o')
%     text(obj.PProcessing.Flux.abs{1,17}(12,3)*1e3,ppval(bh.Extrap_Spline_B,obj.PProcessing.Flux.abs{1,17}(12,3)*1e3)-5.500778160865974e+04,sprintf('%d',ii))
%     hold off

    [MU, KK, ne]  = obj.PProcessing.get_next_permeability(...
        obj.Problem.Materials,obj.Geometry.GeometryElement);
    
    %------------------------------------------------------
    % Update the permeability
    %------------------------------------------------------
    new_material.Permeability = MU;
    obj.Problem = obj.Problem.updateMaterials(new_material);
    [MagMatrix, JacMatrix] = obj.Problem.global_matrix_and_jacobian(KK,obj.PProcessing.Potential);
%     MagMatrix = obj.Problem.global_matrix_contour(KK);
%     JacMatrix = obj.Problem.global_matrix_jacobian(KK,obj.PProcessing.Potential);
    if isempty(obj.Fourier)
        Y_mag = MagMatrix*dataPostProcessing.PHI(SEM_index);
    else
        Y_mag = [MagMatrix*dataPostProcessing.PHI(SEM_index); ...
            zeros(size(obj.Fourier_matrix.Efrequency,1),1)];
    end
    
    obj.Problem = obj.Problem.global_matrix;
    
    s_time = toc;
    
    disp('Concatenating SEM and Fourier Matrices')
    %         E = [obj.Problem.Global_Matrix freq2space; space2freq];
    
    
    
    %         E1 = cat(2,obj.Problem.Global_Matrix,freq2space);
    %         E = cat(1,E1,space2freq);
    
    if isempty(obj.Fourier)
        E = obj.Problem.Global_Matrix;
        JacFin = JacMatrix;
        Y = obj.Problem.Y.vector;
    else
        E(SEM_index,SEM_index)      = obj.Problem.Global_Matrix;
        JacFin(SEM_index,SEM_index) = JacMatrix;
        if ii == 1
            E(end-spr+1:end,end-spc+1:end) = space2freq;
            E(SEM_index(end)-frr+1:SEM_index(end),end-frc+1:end) = freq2space;
            
            JacFin(end-spr+1:end,end-spc+1:end) = space2freq;
            JacFin(SEM_index(end)-frr+1:SEM_index(end),end-frc+1:end) = freq2space;
        end
        Y = cat(2,obj.Problem.Y.vector, zeros(1,size(space2freq,1)));
    end
    fprintf('Computation time for concatenation is %.4f \n',toc - s_time)
    
    
    %     index = [SEM_index fourier_index ];
    %     PHI = zeros(1,(tot_unknowns+numel(fourier_index)));
    
    fprintf('Start solving the linear system at time %.4f \n'...
        ,toc);
    
    
    s_time = toc;
    
    S = sparse(E); JacFin = sparse(JacFin);
%     JacFin = E;
%     JacFin(SEM_index,SEM_index) = JacMatrix; JacFin = sparse(JacFin);
% %     tmp(SEM_index,SEM_index) = JacMatrix - MagMatrix;
%     
%     if ii == 1; dS = S;  prevS = S; else
%         dS = S - prevS;
%     end
%     dS = S; dS(SEM_index,SEM_index) = S(SEM_index,SEM_index) - sparse(MagMatrix);
%     PHI = S\(Y' + Y_mag);% + S\Y_mag;
%     if ii == 1; phi = [PHI; zeros(size(obj.Fourier_matrix.Efrequency,1),1)]; else, phi = PHI; end
if sum(ne)==0
    PHI = S\Y';
    delta_PHI = zeros(size(Y'));
else
    delta_PHI = JacFin\(-S*PHI + Y');
    PHI = PHI + delta_PHI;
%     PHI = S\(Y' + Y_mag);% + S\Y_mag;
%     delta_PHI = prev_PHI - PHI;
%     prev_PHI = PHI;
end
    save(filename,'PHI');
    fprintf('Linear system solved in  %.4f seconds \n',toc - s_time);
    %% pause
    
        dataPostProcessing.PHI          = PHI;
%         dataPostProcessing.PHI(fourier_index) = [];
    
    disp('Evaluating the convergece')
    
    
        obj.PProcessing = PostProcessing(dataPostProcessing);
    
    
    
    
%                 for  k = 1:Nel
%                     err_el(k) = sqrt(sum(((init_phi{k}(:) - obj.PProcessing.Potential{k}(:)).^2).*obj.Geometry.metrics.J{k}(:).*obj.Geometry.metrics.W{k}(:)))./sum(obj.Geometry.metrics.J{k}(:).*obj.Geometry.metrics.W{k}(:));
%                     obj.NonlinearSolver.MU{ii,k} = MU{k};
%                     obj.NonlinearSolver.K{ii,k}  = KK{k};
%                 end
    %                 prev_PHI = PHI;
    err(ii)  = norm(delta_PHI(SEM_index),Inf)/max(abs(PHI(SEM_index)));%/mean(abs(PHI(SEM_index)));
%     prev_PHI = PHI(SEM_index);
%         init_phi = obj.PProcessing.Potential;
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
        ylim([obj.NonlinearSolver.tol*.5 5])
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
    PHI_rem = S(SEM_index,SEM_index)\Y_mag(SEM_index);
else
    PHI_rem = PHI(SEM_index)*0;
    dataPostProcessing.PHI = PHI(SEM_index);
    obj.PProcessing = PostProcessing(dataPostProcessing);
end
dataPostProcessing.PHI = PHI_rem(SEM_index);

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

