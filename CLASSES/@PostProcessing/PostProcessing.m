%==========================================================================
% PostProcessing.m
% Created: 27.10.2017 - 14:49:09
% By: M. Curti
%
% Post processing of the obtained solution
%==========================================================================
classdef PostProcessing < matlab.mixin.SetGet
    %PostProcessing Is taking the solution and extracting and reshaping the
    %solution
    %   Detailed explanation goes here
    
    properties
        PHI
        ProblemData
        xmlContent
        Potential; remPotential; Flux; Permeability
        mappings; metrics; lines; xi; elements; points
        FEM
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = PostProcessing(G, Pr, PHI)
            obj.PHI          = PHI;
            obj.Permeability = Pr.Materials.Permeability;
            obj.ProblemData  = Pr.ProblemData;
            obj.xmlContent   = Pr.xmlContent;
            obj.mappings     = Pr.mappings;
            obj.metrics      = Pr.metrics;
            obj.lines        = G.lines;
            obj.points       = G.points;
            obj.xi           = G.xi;
            obj.elements     = G.elements;
            i_node     = obj.ProblemData.inElementVectorLocation;
            e_node     = obj.ProblemData.lineVectorLocation;
            c_node     = obj.ProblemData.pointVectorLocation;
            
            Nel = obj.ProblemData.Nel;
            obj.Potential   = cell(1,Nel);
            
            for  k = 1:Nel
                % Initiating data
                
                elementsData     = obj.ProblemData.elements;
                
                % Connectivity of elements
                
                if any(k==elementsData.periodic)
                    line_element  = elementsData.periodic_lines(k,:);
                    point_element = elementsData.periodic_points(k,:);
                else
                    line_element  = elementsData.lines(k,:);
                    point_element = elementsData.points(k,:);
                end
               
               
                Potential = zeros(obj.ProblemData.ElementSize{k});
                tmp        = int_el(Potential,'address');
                yel        = tmp{1};
                yb         = tmp{2};
                yc         = tmp{3};
                % Collecting the solutions in each element
                Potential(yel) = obj.PHI(i_node{k});
                
                for ii = 1:4
                    if line_element(ii) < 0
                        Potential(yb{ii}) = ...
                            flipud(obj.PHI(e_node{abs(line_element(ii))}));
                    else
                        Potential(yb{ii}) = ...
                            obj.PHI(e_node{abs(line_element(ii))});
                    end
                end
                
                Potential(yc) = obj.PHI(c_node(point_element));
                obj.Potential{k} = Potential;
                
            end
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Compute the magnetic flux density
        %------------------------------------------------------------------
        
        function obj = compute_B(obj)
            Nel = obj.ProblemData.Nel;
            
            
            obj.Flux.x_comp = cell(1,Nel);
            obj.Flux.y_comp = cell(1,Nel);
            obj.Flux.abs    = cell(1,Nel);
            
            
            % Obtainig the raw results from the B
            for k = 1:Nel
                          
                [Azx, Azy] = ...
                    sem_gradient(obj.Potential{k}, ...
                    obj.metrics.Lx{k}, obj.metrics.Ly{k}, ...
                    obj.metrics.Xcsi{k}, obj.metrics.Xeta{k}, ...
                    obj.metrics.Ycsi{k}, obj.metrics.Yeta{k}, ...
                    obj.metrics.J{k}, obj.ProblemData.ElementSize{k});
                
                obj.Flux.x_comp{k} =  Azy;
                obj.Flux.y_comp{k} = -Azx;
                
                obj.Flux.abs{k} = abs(obj.Flux.x_comp{k} + ...
                                   1i*obj.Flux.y_comp{k});
            end
            
        end
        
         %------------------------------------------------------------------
        % Compute the magnetic flux density
        %------------------------------------------------------------------
        
        function obj = compute_B_axi(obj)
            Nel = obj.ProblemData.Nel;
            
            
            obj.Flux.x_comp = cell(1,Nel);
            obj.Flux.y_comp = cell(1,Nel);
            obj.Flux.abs    = cell(1,Nel);
            
            
            % Obtainig the raw results from the B
            for k = 1:Nel
                          
                [Azx, Azy] = ...
                    sem_gradient(obj.Potential{k}, ...
                    obj.metrics.Lx{k}, obj.metrics.Ly{k}, ...
                    obj.metrics.Xcsi{k}, obj.metrics.Xeta{k}, ...
                    obj.metrics.Ycsi{k}, obj.metrics.Yeta{k}, ...
                    obj.metrics.J{k}, obj.ProblemData.ElementSize{k});
                
                obj.Flux.x_comp{k} =  Azy./obj.mappings.Xm{k};
                obj.Flux.y_comp{k} = -Azx./obj.mappings.Xm{k};
                
                obj.Flux.abs{k} = abs(obj.Flux.x_comp{k} + ...
                                   1i*obj.Flux.y_comp{k});
            end
            
        end
        %------------------------------------------------------------------
        % Plot the solutions on LGL nodes with surf
        %------------------------------------------------------------------
        function plot_surf(obj,varargin)
            Nel = obj.ProblemData.Nel;
                        
            
            if nargin > 1
                mode = varargin{nargin-1};
            else
                mode = 'default';
            end
            
            switch mode
                case 'default'
                    for k = 1:Nel
                        surf(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                            obj.Flux.abs{k},'edgecolor','none');
                    end
                case 'potential'
                    for k = 1:Nel
                        surf(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                            obj.Potential{k},'edgecolor','none');
                    end
            end
        end
        %==================================================================
        
        %------------------------------------------------------------------
        % Plot the solution on the LGL nodes with contour
        %------------------------------------------------------------------
        function plot_contour(obj,varargin)
            Nel = obj.ProblemData.Nel;
            
            % Computing the flux lines
            if nargin == 2
                flux_lines = ...
                            linspace(min(varargin{1}),max(varargin{1}),20);
            else
                flux_lines = linspace(min(obj.PHI(:)),max(obj.PHI(:)),20);
            end
            
            % Plot contours
            for k = 1:Nel
                contour(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                    obj.Potential{k},flux_lines,'linewidth',1);
            end
        end
        
        %------------------------------------------------------------------
        % Plot the solution on the LGL nodes with contour
        %------------------------------------------------------------------
        function plot_contour_axi(obj,varargin)
            Nel = obj.ProblemData.Nel;
            
            % Computing the flux lines
            if nargin == 2
                flux_lines = ...
                            linspace(min(varargin{1}),max(varargin{1}),20);
            else
                flux_lines = linspace(min(obj.PHI(:)),max(obj.PHI(:)),20);
            end
            
            % Plot contours
            for k = 1:Nel
                contour(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                    obj.Potential{k}*1./obj.mappings.Xm{k},flux_lines);
            end
        end
        
        %------------------------------------------------------------------
        % Plot a quantity on SEM nodes
        %------------------------------------------------------------------
        function plot_surf_var(obj,var,Nel)
            if nargin ==2
            Nel = 1:obj.ProblemData.Nel;
            end
            
                       
                    for k = Nel
                        surf(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                            var{1,k},'edgecolor','none');
                    end
        end
        
        %------------------------------------------------------------------
        % Plot the a quantity on the LGL nodes with contour
        %------------------------------------------------------------------
        function plot_contour_var(obj,var)
            Nel = obj.ProblemData.Nel;
            
            % findig the min and max
            min_var = 0; max_var = 0;
            for k = 1:Nel
                min_var = min([min_var min(var{k}(:))]);
                max_var = max([max_var max(var{k}(:))]);
            end
            
            flux_lines = linspace(min_var, max_var, 20);
            
            % Plot contours
            for k = 1:Nel
            contour(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                    var{k},flux_lines);
               
            end
        end
        
        %------------------------------------------------------------------
        % Plot the solution on FEM nodes
        %------------------------------------------------------------------
        
        function plot_trisurf(obj,varargin)
            Nel      = obj.ProblemData.Nel;
            dataFlux = varargin{1};
               for k = 1:Nel 
                % initialisation of data
                x_mesh  = dataFlux.x_mesh{k};
                y_mesh  = dataFlux.y_mesh{k};
                Xm      = obj.mappings.Xm{k};
                Ym      = obj.mappings.Ym{k};
                TRI     = dataFlux.TRI{k};
                PHI_FEM = dataFlux.values{k};
                
                    % Interpolation matrix for the k th element
                   lxy        = PolynomialInterp12D(x_mesh,y_mesh,Xm...
                                                              *1e3,Ym*1e3);
                   phi_SEM    = lxy*obj.Potential{k}(:);
                   trisurf(TRI,x_mesh*1e-3,y_mesh*1e-3,...
                   abs(phi_SEM-PHI_FEM),'edgecolor','none')
               end
        end
        %==================================================================
        %------------------------------------------------------------------
        % Interpolate to all FEM nodes
        %------------------------------------------------------------------

        SEM = InterpolateOnFEMNodes(obj,allFEM)
        %------------------------------------------------------------------
        % Function which combines all the points from the elements in one
        % set of vectors
        function allSEM = allSEMinOne(obj,varargin)
            % initialisation of the variables
            Nel              = obj.ProblemData.Nel;
            allSEM.x         = zeros(obj.ProblemData.GridPoints,1);
            allSEM.y         = zeros(obj.ProblemData.GridPoints,1);
            allSEM.potential = zeros(obj.ProblemData.GridPoints,1);
            allSEM.FluxAbs   = zeros(obj.ProblemData.GridPoints,1);
            
            % Checking if there are some specified elements
            if nargin==2
                this_elements = varargin{1};
            else
                this_elements = 1:Nel;
            end
            % main for loop
            this_index = 0;
            for k = this_elements
                ElementSize = prod(obj.ProblemData.ElementSize{k});
                allSEM.x((1:ElementSize) + this_index) = ...
                                                     obj.mappings.Xm{k}(:);
                allSEM.y((1:ElementSize) + this_index) = ...
                                                     obj.mappings.Ym{k}(:);
                allSEM.potential((1:ElementSize) + this_index) = ...
                                                     obj.Potential{k}(:);
                allSEM.FluxAbs((1:ElementSize) + this_index) = ...
                                                     obj.Flux.abs{k}(:);
                this_index = this_index + ElementSize;
            end
            
            % Excluding the identical points
            [~,IA,~] = unique(round(allSEM.x + ...
                           1i*allSEM.y,12),'stable');
            allSEM.x         = allSEM.x(IA);
            allSEM.y         = allSEM.y(IA);
            allSEM.potential = allSEM.potential(IA);
            allSEM.FluxAbs   = allSEM.FluxAbs(IA);
        end
        %------------------------------------------------------------------
        % Function that returns the updated permeability according to the
        % given BH rule and the source as a result of remanence
        [Permeability, K, nonlin_elem ]= ...
                            get_next_permeability(obj,Materials,xmlContent)
        %------------------------------------------------------------------
        % Funciton to solve the flux linkage
        %------------------------------------------------------------------
         [FL, remFL] = flux_linkage(obj,elements)
        
        %------------------------------------------------------------------
        % Function to solve the force on the elements
        %------------------------------------------------------------------
        function [Fxt, Fyt] = element_force(obj, elements)
            element_lines = obj.elements.lines(elements,:);
            involved_lines = unique(element_lines);
            boundary_lines = zeros(size(involved_lines));
            mu0 = pi*4e-7;
            k = 1;
            
            % Selecting the boundary lines
            for ll = involved_lines'
                num_el = numel(find(element_lines==ll));
                if num_el == 1
                    boundary_lines(k) = ll;
                    k = k + 1;
                end
                boundary_lines(boundary_lines == 0) = [];
            end
            
            % Finding the outside elements
            outside_elements = zeros(numel(boundary_lines),2);
            inside_elements  = zeros(numel(boundary_lines),2);
            shl = obj.lines.shared_line(boundary_lines,:);
            
            Fx = zeros(1,numel(boundary_lines));
            Fy = zeros(1,numel(boundary_lines));
            
            for ll = 1:length(outside_elements)
                % Filling in the list of outside elements
                if any(shl(ll,1)==elements)
                    outside_elements(ll,1) = shl(ll,2);
                    outside_elements(ll,2) = shl(ll,4);
                    inside_elements(ll,1)  = shl(ll,1);
                    inside_elements(ll,2)  = shl(ll,3);
                else
                    outside_elements(ll,1) = shl(ll,1);
                    outside_elements(ll,2) = shl(ll,3);
                    inside_elements(ll,1)  = shl(ll,2);
                    inside_elements(ll,2)  = shl(ll,4);
                end
                % Checking if the element is outside of the domain
                if outside_elements(ll,1) == 0
                    error('The element must not touch the boundary')
                end
                
                % Extracting the physical and geometrical parameters
                % The quadratures
                l_id_out = outside_elements(ll,2);
                e_id_out = outside_elements(ll,1);
                l_id_in = inside_elements(ll,2);
                e_id_in = inside_elements(ll,1);
                lg_id = boundary_lines(ll);
                if l_id_out==1 || l_id_out==3
                    neu_tmp  = obj.metrics.Neu2{e_id_out};
                    if l_id_out == 1
                        n_tmpx   = -obj.metrics.N2{1,e_id_in};
                        n_tmpy   = -obj.metrics.N2{2,e_id_in};
                    else
                        n_tmpx   =  obj.metrics.N2{1,e_id_in};
                        n_tmpy   =  obj.metrics.N2{2,e_id_in};
                    end
                else
                    neu_tmp  = obj.metrics.Neu1{e_id_out};
                    if l_id_out == 2
                        n_tmpx   =  obj.metrics.N1{1,e_id_in};
                        n_tmpy   =  obj.metrics.N1{2,e_id_in};
                    else
                        n_tmpx   =  -obj.metrics.N1{1,e_id_in};
                        n_tmpy   =  -obj.metrics.N1{2,e_id_in};
                    end
                end
                tmp    = int_el(neu_tmp,'lines_address');
                l_pts_out = tmp{l_id_out};
                l_pts_in = tmp{l_id_in};
                
                wl = neu_tmp(l_pts_out).*obj.xi.w_for_all_lines{lg_id}';
                nx = n_tmpx(l_pts_in);
                ny = n_tmpy(l_pts_in);
                Bx = obj.Flux.x_comp{e_id_out}(l_pts_out);
                By = obj.Flux.y_comp{e_id_out}(l_pts_out);
                B  = obj.Flux.abs{e_id_out}(l_pts_out);
                Fx(ll) = ((Bx.^2-B.^2*.5).*nx + Bx.*By.*ny)*wl'/mu0;
                                       
                Fy(ll) = ((By.^2-B.^2*.5).*ny + Bx.*By.*nx)*wl'/mu0;
            end
            Fxt = sum(Fx); Fyt = sum(Fy);
        end
    end
    
end



function [fx, fy] = sem_gradient(f, Lx, Ly, Xxi, Xeta, Yxi, Yeta, J, ...
                                                               ElementSize)
%    fxi  = reshape(Lx*f(:),ElementSize(1),ElementSize(2));
%                    
%    feta = reshape(Ly*f(:),ElementSize(1),ElementSize(2));
%                    
%                    
%    fx = (Yeta.*fxi - Yxi.*feta)./J;
%    
%    fy = (-Xeta.*fxi + Xxi.*feta)./J;
% The same thing in the matrix form
   fx = reshape(diag(1./J(:))*(diag(Yeta(:))*Lx - diag(Yxi(:))*Ly)*f(:),ElementSize(1),ElementSize(2));
   
   fy = reshape(diag(1./J(:))*(-diag(Xeta(:))*Lx + diag(Xxi(:))*Ly)*f(:),ElementSize(1),ElementSize(2));
end