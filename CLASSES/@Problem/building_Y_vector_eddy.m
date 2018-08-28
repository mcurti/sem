%------------------------------------------------------------------
% Building the right hand side vector
%------------------------------------------------------------------

function vector = building_Y_vector_eddy(obj)
Nel    = obj.ProblemData.Nel;

index_Elements = obj.ProblemData.inElementVectorLocation;
index_lines    = obj.ProblemData.lineVectorLocation;
index_points   = obj.ProblemData.pointVectorLocation;
elementsData   = obj.ProblemData.elements;

vector = ...
    zeros(1,obj.ProblemData.pointVectorLocation(end));
for k = 1:Nel
    EL = obj.Y.Conductivity{k};
    % Getting the connectivity of the lines and points of
    % k-th element
    if any(k==elementsData.periodic)
        line_element  = elementsData.periodic_lines(k,:);
        point_element = elementsData.periodic_points(k,:);
        line_element_no_abs = line_element;
    else
        line_element  = abs(elementsData.lines(k,:));
        line_element_no_abs = elementsData.lines(k,:);
        point_element = elementsData.points(k,:);
    end
    
    % Boundary Conditions
    dir_points = eval(['[',xml2matlab(obj.xmlContent...
        ,'BoundaryConditions',0,'dir_points','Attribute'),'];']);
    dir_lines  = eval(['[',xml2matlab(obj.xmlContent...
        ,'BoundaryConditions',0,'dir_lines','Attribute'),'];']);
    % Filling the inElement content
    vector(index_Elements{k}) = int_el(EL,'intElement');
    
    % Filling the line content
    for ii = 1:4
        %if any(line_element(ii)==dir_lines)
        %else
            tmp_line = ...
                cell2mat(int_el(EL,ii,'lines'));
            if line_element_no_abs(ii)<0
                vector(index_lines{line_element(ii)})=...
                    vector(index_lines{line_element(ii)})+...
                    fliplr(tmp_line);
                
            else
                vector(index_lines{line_element(ii)})=...
                    vector(index_lines{line_element(ii)})+...
                    tmp_line;
            end
        %end
    end
    
    % Filling the points content
    for ii = 1:4
        %if any(point_element(ii)==dir_points)
        %else
            vector(index_points(point_element(ii))) = ...
                vector(index_points(point_element(ii))) + ...
                int_el(EL,ii,'points');
        %end
    end
end

% Filling the imposed Flux boundary condition

end
%==================================================================
