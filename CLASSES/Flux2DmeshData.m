classdef Flux2DmeshData
    %Flux2DmeshData extracting and preparing for visualisation or post
    %processing of data exported from Flux2D 
    %   Detailed explanation goes here
    
    properties
        filenames
        x_mesh; y_mesh; values_mesh; TRI_mesh;
        total_nodes
        
    end
    
    methods
        %Constructor
        function obj = Flux2DmeshData(filenames)
            
            if iscell(filenames)==false
                filenames = {filenames};
            end
            obj.filenames = filenames;
            files_number = numel(filenames);
            
            obj.x_mesh      = cell(1,files_number);
            obj.y_mesh      = cell(1,files_number);
            obj.values_mesh = cell(1,files_number);
            obj.TRI_mesh    = cell(1,files_number);
            obj.total_nodes = 0;
            
            for k = 1:files_number
                raw_data           = import_spatial_quantity(filenames{k});
                
                obj.x_mesh{k}      = raw_data(:,1);
                obj.y_mesh{k}      = raw_data(:,2);
                obj.values_mesh{k} = raw_data(:,3);
                obj.TRI_mesh{k}    = delaunay(obj.x_mesh{k},obj.y_mesh{k});
                obj.total_nodes    = obj.total_nodes + ...
                                                      numel(obj.x_mesh{k});
            end
        end
        %==================================================================
        
                
        %------------------------------------------------------------------
        % FEMfilesInOne is grouping the points and values from all files
        % from the class into one vector. It also excludes the triangles
        % which look odd.
        %------------------------------------------------------------------
        
        function allFEM = FEMfilesInOne(obj)
            
            % Pile all the files data in one vector
            filenumber = numel(obj.filenames);
            
            allFEM.all_x_nodes = zeros(obj.total_nodes,1);
            allFEM.all_y_nodes = zeros(obj.total_nodes,1);
            allFEM.all_values  = zeros(obj.total_nodes,1);
            nodes_counter       = 0;
            
            for k = 1:filenumber
                file_nodes = numel(obj.x_mesh{k});
                
                allFEM.all_x_nodes((1:file_nodes) + nodes_counter) =...
                    obj.x_mesh{k};
                
                allFEM.all_y_nodes((1:file_nodes) + nodes_counter) =...
                    obj.y_mesh{k};
                
                allFEM.all_values((1:file_nodes) + nodes_counter)  =...
                    obj.values_mesh{k};
                
                nodes_counter = nodes_counter + file_nodes;
            end
            
            % Excluding the similar points
            
            [~,IA,~] = unique(allFEM.all_x_nodes + ...
                           1i*allFEM.all_y_nodes,'stable');
            allFEM.all_x_nodes = allFEM.all_x_nodes(IA);
            allFEM.all_y_nodes = allFEM.all_y_nodes(IA);
            allFEM.all_values  = allFEM.all_values(IA);
            
            % Building the mesh
            allFEM.TRI = delaunay(allFEM.all_x_nodes,allFEM.all_y_nodes);
            % Measuring the sides of the triangles
            sides = zeros(size(allFEM.TRI));
            x = allFEM.all_x_nodes; y = allFEM.all_y_nodes;
            tri = allFEM.TRI;
            
            sides(:,1) = abs(x(tri(:,1)) + 1i*y(tri(:,1)) - ...
                            (x(tri(:,2)) + 1i*y(tri(:,2))));
            sides(:,2) = abs(x(tri(:,2)) + 1i*y(tri(:,2)) - ...
                            (x(tri(:,3)) + 1i*y(tri(:,3))));
            sides(:,3) = abs(x(tri(:,3)) + 1i*y(tri(:,3)) - ...
                            (x(tri(:,1)) + 1i*y(tri(:,1))));
            sides_ratio = max(sides,[],2)./min(sides,[],2);
            
            % Measuring the area of the triangles
            AR      = pdetrg([x y]',tri')';
            meanAR  = mean(AR);
            
            ratio_criterea = sides_ratio > 2;
            area_criterea  = AR > meanAR*2;
            
            sum_criterea = bitor(area_criterea,ratio_criterea);
            
            allFEM.TRI(sum_criterea,:) = [];
        end
        
        % Plot mesh function
        function plot_mesh(obj,varargin)
            files_number = numel(obj.filenames);
            
            mode = varargin{nargin-1};
            switch mode
                %========================
                case 'all'
                %========================
                    % In the case the data from all mesh elements has to be
                    % plotted
                    for k = 1:files_number
                        hold on
                        % Plot the element
                        triplot(obj.TRI_mesh{k},obj.x_mesh{k},...
                            obj.y_mesh{k})
                        hold off
                    end
                    % In the case that non square elements are plotted,
                    % then it is convenient to put a white rectangle in the
                    % background of some elements, so that the nasty lines
                    % from the unwanted triangles are not visible
                %=========================
                case 'all_wite_rectangles'
                %=========================
                    % The vector which indicates which element will have a
                    % white background.
                    rectangles = varargin{1};
                    for k = 1:files_number
                        hold on
                        % Plot an white rectangle in the back of each mesh
                        % element
                        if any(k == rectangles)
                            x_rec_bot = min(obj.x_mesh{k}); 
                            x_rec_top = max(obj.x_mesh{k});
                            
                            y_rec_bot = min(obj.y_mesh{k}); 
                            y_rec_top = max(obj.y_mesh{k});
                            
                            rectangle('position',[x_rec_bot y_rec_bot ...
                            x_rec_top - x_rec_bot y_rec_top - y_rec_bot...
                            ],'facecolor','w')
                        end
                        % Plot the element
                        triplot(obj.TRI_mesh{k},obj.x_mesh{k},...
                            obj.y_mesh{k},'color','k')
                        hold off
                    end
                    
            end
        end
        
        %===============================================================
        % add your own function
    end
    
end

