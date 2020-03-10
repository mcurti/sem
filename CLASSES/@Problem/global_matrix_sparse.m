function obj = global_matrix_sparse(obj)
% Data initialisation
Nel      = obj.ProblemData.Nel;
Nl       = obj.ProblemData.Nl;
Np       = obj.ProblemData.Np;
unknowns = obj.ProblemData.pointVectorLocation(end);

% E_sum  = zeros(obj.ProblemData.pointVectorLocation(end));
% I      = eye(obj.ProblemData.pointVectorLocation(end));

% Initalisation of the sparse matrix

nz = round(unknowns^2/Nel);

rows = zeros(nz,1);
cols = zeros(nz,1);
v  = zeros(nz,1);
current_index = 0;
% Boundary Conditions
dir_points = eval(['[',xml2matlab(obj.xmlContent...
    ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
    ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
per_points = eval(['[',xml2matlab(obj.xmlContent...
    ,'BoundaryConditions',0,'periodic_points','Attribute'),'];']);
per_lines  = eval(['[',xml2matlab(obj.xmlContent...
    ,'BoundaryConditions',0,'periodic_lines','Attribute'),'];']);


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
    E = element_matrix(dataElement);
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
    end
    
    %------------------------------------------------------------
    %E_sum(i_node,i_node) =  E(yel,yel);
    [tmp_rows, tmp_cols] = meshgrid(i_node,i_node); tmp_v = E(yel,yel);
    
    tmp_length = numel(tmp_rows);
    index_vector = (1:tmp_length) + current_index;
    
    rows(index_vector) = tmp_rows(:);
    cols(index_vector) = tmp_cols(:);
    v   (index_vector) = tmp_v(:);
    
    current_index = current_index + tmp_length;
    %=================================================================
    jj = abs(line_element); ll = point_element;
    % inElement Points
    
    %E_sum(i_node,c_node(ll)) = E_sum(i_node,c_node(ll)) + E(yel,yc);
    %==================================================================
    [tmp_rows, tmp_cols] = meshgrid(i_node,c_node(ll)); tmp_v = E(yel,yc);
    
    tmp_length = numel(tmp_rows);
    index_vector = (1:tmp_length) + current_index;
    
    rows(index_vector) = tmp_rows(:);
    cols(index_vector) = tmp_cols(:);
    v   (index_vector) = tmp_v(:);
    
    current_index = current_index + tmp_length;
    %=================================================================
    
    % Points inElement
    %E_sum(c_node(ll),i_node) = E_sum(c_node(ll),i_node) + E(yc,yel);
    %==================================================================
    [tmp_rows, tmp_cols] = meshgrid(c_node(ll),i_node); tmp_v = E(yc,yel);
    
    tmp_length = numel(tmp_rows);
    index_vector = (1:tmp_length) + current_index;
    
    rows(index_vector) = tmp_rows(:);
    cols(index_vector) = tmp_cols(:);
    v   (index_vector) = tmp_v(:);
    
    current_index = current_index + tmp_length;
    %=================================================================
    
    % Points Points
    % E_sum(c_node(ll),c_node(ll)) = ...
    % E_sum(c_node(ll),c_node(ll)) + E(yc,yc);
    %==================================================================
    [tmp_rows, tmp_cols] = meshgrid(c_node(ll),c_node(ll)); tmp_v = E(yc,yc);
    
    tmp_length = numel(tmp_rows);
    index_vector = (1:tmp_length) + current_index;
    
    rows(index_vector) = tmp_rows(:);
    cols(index_vector) = tmp_cols(:);
    v   (index_vector) = tmp_v(:);
    
    current_index = current_index + tmp_length;
    %=================================================================
    
    
    
    for ii = 1:4
        % element edges
        % E_sum(i_node,e_node{jj(ii)}) =   ...
        % E_sum(i_node,e_node{jj(ii)}) + E(yel,yb{ii});
        %==================================================================
        [tmp_rows, tmp_cols] = meshgrid(i_node,e_node{jj(ii)}); tmp_v = E(yel,yb{ii});
        
        tmp_length = numel(tmp_rows);
        index_vector = (1:tmp_length) + current_index;
        
        rows(index_vector) = tmp_rows(:);
        cols(index_vector) = tmp_cols(:);
        v   (index_vector) = tmp_v(:);
        
        current_index = current_index + tmp_length;
        %=================================================================
        
        % edges element
        
        % E_sum(e_node{jj(ii)},i_node) =   ...
        % E_sum(e_node{jj(ii)},i_node) + E(yb{ii},yel);
        %==================================================================
        [tmp_rows, tmp_cols] = meshgrid(e_node{jj(ii)},i_node); tmp_v = E(yb{ii},yel);
        
        tmp_length = numel(tmp_rows);
        index_vector = (1:tmp_length) + current_index;
        
        rows(index_vector) = tmp_rows(:);
        cols(index_vector) = tmp_cols(:);
        v   (index_vector) = tmp_v(:);
        
        current_index = current_index + tmp_length;
        %=================================================================
        
        % edges edges
        for n = 1:4
            % E_sum(e_node{jj(ii)},e_node{jj(n)}) = E_sum(...
            % e_node{jj(ii)},e_node{jj(n)}) + E(yb{ii},yb{n});
            %==================================================================
            [tmp_rows, tmp_cols] = meshgrid(e_node{jj(ii)},e_node{jj(n)}); tmp_v = E(yb{ii},yb{n});
            
            tmp_length = numel(tmp_rows);
            index_vector = (1:tmp_length) + current_index;
            
            rows(index_vector) = tmp_rows(:);
            cols(index_vector) = tmp_cols(:);
            v   (index_vector) = tmp_v(:);
            
            current_index = current_index + tmp_length;
            %=================================================================
        end
        
        % edges corners
        % E_sum(e_node{jj(ii)},c_node(ll)) = ...
        % E_sum(e_node{jj(ii)},c_node(ll)) + E(yb{ii},yc);
        %==================================================================
        [tmp_rows, tmp_cols] = meshgrid(e_node{jj(ii)},c_node(ll)); tmp_v = E(yb{ii},yc);
        
        tmp_length = numel(tmp_rows);
        index_vector = (1:tmp_length) + current_index;
        
        rows(index_vector) = tmp_rows(:);
        cols(index_vector) = tmp_cols(:);
        v   (index_vector) = tmp_v(:);
        
        current_index = current_index + tmp_length;
        %=================================================================
        
        % corners edges
        % E_sum(c_node(ll),e_node{jj(ii)}) = ...
        % E_sum(c_node(ll),e_node{jj(ii)}) + E(yc,yb{ii});
        %==================================================================
        [tmp_rows, tmp_cols] = meshgrid(c_node(ll),e_node{jj(ii)}); tmp_v = E(yc,yb{ii});
        
        tmp_length = numel(tmp_rows);
        index_vector = (1:tmp_length) + current_index;
        
        rows(index_vector) = tmp_rows(:);
        cols(index_vector) = tmp_cols(:);
        v   (index_vector) = tmp_v(:);
        
        current_index = current_index + tmp_length;
        %=================================================================
    end
    
    
    
end


% dirichlet
% points
for k = 1:Np
    if any(k==dir_points)
%         E_sum(c_node(k),:)  = I(c_node(k),:);
        v(rows==c_node(k))  = 0;
        v(rows==c_node(k)& cols==c_node(k))  = 1;
    end
end

% lines
for k = 1:Nl
    if any(k==dir_lines)
        %         E_sum(e_node{k},:) = I(e_node{k},:);
        for kk = 1:numel(e_node{k})
            v(rows==e_node{k}(kk))  = 0;
            v(rows==e_node{k}(kk)& cols==e_node{k}(kk))  = 1;
        end
    end
end

% periodic
% points
for k = 1:Np
    if any(k==per_points)
        %                    E_sum(:,c_node(k))  = 0;
%         E_sum(c_node(k),:)  = I(c_node(k),:);
        v(rows==c_node(k))  = 0;
        v(rows==c_node(k)& cols==c_node(k))  = 1;
    end
end

% lines
for k = 1:Nl
    if any(k==per_lines)
        %                    E_sum(:,e_node{k}) = 0;
%         E_sum(e_node{k},:) = I(e_node{k},:);
        for kk = 1:numel(e_node{k})
            v(rows==e_node{k}(kk))  = 0;
            v(rows==e_node{k}(kk)& cols==e_node{k}(kk))  = 1;
        end
    end
end
cols(rows==0) = []; v(rows==0) = []; rows(rows==0) = [];
E_sum = sparse(rows,cols,v);
obj.Global_Matrix = E_sum;
end

%==================================================================

