
%------------------------------------------------------------------
% Interpolate to all FEM nodes
%------------------------------------------------------------------

function  SEM = InterpolateOnFEMNodes(obj,allFEM)

% Initialising the data
SEM.potential = zeros(size(allFEM.all_x_nodes));
SEM.abs_Flux  = zeros(size(allFEM.all_x_nodes));
Nel = obj.ProblemData.Nel;
el_list           = 1:Nel;
obj.FEM.x_mesh    = cell(1,Nel);
obj.FEM.y_mesh    = cell(1,Nel);
obj.FEM.TRI       = cell(1,Nel);
obj.FEM.potential = cell(1,Nel);
obj.FEM.Flux_abs  = cell(1,Nel);
s                 = cell(2,4);
options = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-8);

% Element list
try
    el_list(allFEM.noelements) = [];
catch
    disp('interpolating with all elements')
end

% the loop for all FEM nodes
for ii = 1:length(allFEM.all_x_nodes)
    
    x_fem = allFEM.all_x_nodes(ii);
    y_fem = allFEM.all_y_nodes(ii);
    interpolated = false;
    
    
    % Check if init values are present
    if exist('csi_init','var')
    else
        csi_init = 1-2*rand(1); eta_init = 1-2*rand(1);
    end
    % the loop for each element
    try
        tmp = allFEM.int;
        el  = tmp.element(ii);
        csi_i = tmp.csi(ii);
        eta_i = tmp.eta(ii);
        connectivity = abs(obj.elements.lines(el,:));
        csi = obj.xi.csi_for_all_lines{connectivity(1)};
        eta = obj.xi.csi_for_all_lines{connectivity(2)};
        [X, Y] = meshgrid(csi, eta);
        lxy        = PolynomialInterp12D(csi_i,eta_i,X,Y);
        SEM.potential(ii)   = lxy*obj.Potential{el}(:);
        SEM.abs_Flux(ii)    = lxy*obj.Flux.abs{el}(:)*1e3;
    catch
        for jj = el_list
            clc
            fprintf('interpolation on node %d from %d in the element %d\n',ii,...
                length(allFEM.all_x_nodes),jj);
            connectivity = abs(obj.elements.lines(jj,:));
            
            for k = 1:4
                s{1,k} = obj.lines.vector{connectivity(k)}(1,:);
                s{2,k} = obj.lines.vector{connectivity(k)}(2,:);
            end
            
            % Check if the point is close to the element
            
            max_dx = max(max(obj.mappings.Xm{jj}(:)) - ...
                min(obj.mappings.Xm{jj}(:)));
            
            max_dy = max(max(obj.mappings.Ym{jj}(:)) - ...
                min(obj.mappings.Ym{jj}(:)));
            
            min_dx_FEM = min(abs(x_fem - obj.mappings.Xm{jj}(:)));
            min_dy_FEM = min(abs(y_fem - obj.mappings.Ym{jj}(:)));
            
            
            if max_dx > min_dx_FEM && max_dy > min_dy_FEM
                
                % Axes in the computational domain
                csi = obj.xi.csi_for_all_lines{connectivity(1)};
                eta = obj.xi.csi_for_all_lines{connectivity(2)};
                
                myfun = @(x)  [x_fem; y_fem] - ...
                    InterpolationPoints(csi,eta,x(1:numel(x_fem)),...
                    x(numel(x_fem)+1:end), s);
                
                if ii==1
                    [out,~,exitflag,~] = ...
                        fsolve(myfun,1-2*rand(1,2),options);
                else
                    [out,~,exitflag,~] = ...
                        fsolve(myfun,[csi_init, eta_init],options);
                end
                
                
                
                csi_i = out(1); eta_i = out(2);
                
                if abs(round(csi_i,6)) <= 1 && ...
                        abs(round(eta_i,6)) <= 1 && exitflag > 0
                    [X, Y] = meshgrid(csi, eta);
                    lxy        = PolynomialInterp2D(csi_i,eta_i,X,Y);
                    
                    SEM.potential(ii)   = lxy*obj.Potential{jj}(:);
                    SEM.abs_Flux(ii)    = lxy*obj.Flux.abs{jj}(:);
                    SEM.int.element(ii) = jj;
                    SEM.int.csi(ii)     = csi_i;
                    SEM.int.eta(ii)     = eta_i;
                    interpolated = true;
                    csi_init = out(1); eta_init = out(2);
                end
            end
            
            % Break the element loop if interpolated
            if interpolated
                break
            end
        end
    end
end
end
