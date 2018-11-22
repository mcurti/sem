function [G, Ja] = global_matrix_and_jacobian(obj,kk,Potential)
% Data initialisation
Nel    = obj.ProblemData.Nel;
Nl     = obj.ProblemData.Nl;
Np     = obj.ProblemData.Np;

E_sum  = zeros(obj.ProblemData.pointVectorLocation(end));
J_sum  = zeros(obj.ProblemData.pointVectorLocation(end));

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
dataElement.bh = bh_class(2);

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
    dataElement.a    = Potential{k}(:);
    dataElement.kk = kk{k};
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
    [E, J] = element_matrix_GJ(dataElement);
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
        %---------------------------------------------
        
        LL = J(inversed_line,inversed_line);
        EL = J(:,inversed_line);
        LE = J(inversed_line,:);
        J(inversed_line,:) = flipud(LE);
        J(:,inversed_line) = fliplr(EL);
        J(inversed_line,inversed_line) = rot90(LL,2);
    end
    
    %------------------------------------------------------------
    E_sum(i_node,i_node) =  E(yel,yel);
    J_sum(i_node,i_node) =  J(yel,yel);
    %=================================================================
    jj = abs(line_element); ll = point_element;
    % inElement Points
    
    E_sum(i_node,c_node(ll)) = E_sum(i_node,c_node(ll)) + E(yel,yc);
    J_sum(i_node,c_node(ll)) = J_sum(i_node,c_node(ll)) + J(yel,yc);
    
    % Points inElement
    E_sum(c_node(ll),i_node) = E_sum(c_node(ll),i_node) + E(yc,yel);
    J_sum(c_node(ll),i_node) = J_sum(c_node(ll),i_node) + J(yc,yel);
    
    % Points Points
    E_sum(c_node(ll),c_node(ll)) = ...
        E_sum(c_node(ll),c_node(ll)) + E(yc,yc);
    
    J_sum(c_node(ll),c_node(ll)) = ...
        J_sum(c_node(ll),c_node(ll)) + J(yc,yc);
    
    
    
    for ii = 1:4
        % element edges
        E_sum(i_node,e_node{jj(ii)}) =   ...
            E_sum(i_node,e_node{jj(ii)}) + E(yel,yb{ii});
        J_sum(i_node,e_node{jj(ii)}) =   ...
            J_sum(i_node,e_node{jj(ii)}) + J(yel,yb{ii});
        
        % edges element
        
        E_sum(e_node{jj(ii)},i_node) =   ...
            E_sum(e_node{jj(ii)},i_node) + E(yb{ii},yel);
        J_sum(e_node{jj(ii)},i_node) =   ...
            J_sum(e_node{jj(ii)},i_node) + J(yb{ii},yel);
        
        % edges edges
        for n = 1:4
            E_sum(e_node{jj(ii)},e_node{jj(n)}) = E_sum(...
                e_node{jj(ii)},e_node{jj(n)}) + E(yb{ii},yb{n});
            J_sum(e_node{jj(ii)},e_node{jj(n)}) = J_sum(...
                e_node{jj(ii)},e_node{jj(n)}) + J(yb{ii},yb{n});
        end
        
        % edges corners
        E_sum(e_node{jj(ii)},c_node(ll)) = ...
            E_sum(e_node{jj(ii)},c_node(ll)) + E(yb{ii},yc);
        J_sum(e_node{jj(ii)},c_node(ll)) = ...
            J_sum(e_node{jj(ii)},c_node(ll)) + J(yb{ii},yc);
        
        % corners edges
        E_sum(c_node(ll),e_node{jj(ii)}) = ...
            E_sum(c_node(ll),e_node{jj(ii)}) + E(yc,yb{ii});
        J_sum(c_node(ll),e_node{jj(ii)}) = ...
            J_sum(c_node(ll),e_node{jj(ii)}) + J(yc,yb{ii});
    end
    
end

% dirichlet
% points
for k = 1:Np
    if any(k==dir_points)
        E_sum(c_node(k),:)  = I(c_node(k),:);
        J_sum(c_node(k),:)  = I(c_node(k),:);
    end
end

% lines
for k = 1:Nl
    if any(k==dir_lines)
        E_sum(e_node{k},:) = I(e_node{k},:);
        J_sum(e_node{k},:) = I(e_node{k},:);
    end
end

% periodic
% points
for k = 1:Np
    if any(k==per_points)
        %                    E_sum(:,c_node(k))  = 0;
        E_sum(c_node(k),:)  = I(c_node(k),:);
        J_sum(c_node(k),:)  = I(c_node(k),:);
    end
end

% lines
for k = 1:Nl
    if any(k==per_lines)
        %                    E_sum(:,e_node{k}) = 0;
        E_sum(e_node{k},:) = I(e_node{k},:);
        J_sum(e_node{k},:) = I(e_node{k},:);
    end
end

G = E_sum;
Ja = J_sum;
end