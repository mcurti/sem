%==========================================================================
% PostProcessing.m
% Created: 27.10.2017 - 14:49:09
% By: M. Curti
%
% Post processing of the obtained solution
%==========================================================================
classdef PostProcessing
    %PostProcessing Is taking the solution and extracting and reshaping the
    %solution
    %   Detailed explanation goes here
    
    properties
        PHI
        ProblemData
        xmlContent
        Potential; Flux; Permeability
        mappings; metrics; lines; xi; elements; points
        FEM
    end
    
    methods
        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        function obj = PostProcessing(dataPostProcessing)
            obj.PHI          = dataPostProcessing.PHI;
            obj.Permeability = dataPostProcessing.Permeability;
            obj.ProblemData  = dataPostProcessing.ProblemData;
            obj.xmlContent   = dataPostProcessing.xmlContent;
            obj.mappings     = dataPostProcessing.mappings;
            obj.metrics      = dataPostProcessing.metrics;
            obj.lines        = dataPostProcessing.lines;
            obj.points       = dataPostProcessing.points;
            obj.xi           = dataPostProcessing.xi;
            obj.elements     = dataPostProcessing.elements;
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
%             Nl  = obj.ProblemData.Nl;
%             Np  = obj.ProblemData.Np;
            
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
%             
%             % Post processing the derivatives based on boundary conditions
%             
%             % Adjusting the values on the lines
%             shared_line = abs(obj.lines.shared_line);
%             for k = 1:Nl
%                 % Filter the lines which do not share two elements
%                 if isempty(find(shared_line(k,:)==0,1)) == true
%                     % Get the line values on both elements
%                     el1 = shared_line(k,1); el2 = shared_line(k,2);
%                     sd1 = shared_line(k,3); sd2 = shared_line(k,4);
%                     
%                     tmp1 = int_el(obj.Flux.x_comp{el1},'address');
%                     tmp2 = int_el(obj.Flux.x_comp{el2},'address');
%                     
%                     b1 = tmp1{2}; b2 = tmp2{2};
%                     
%                     % Permeability values
%                     line1_mu = obj.Permeability{el1}(b1{sd1});
%                     line2_mu = obj.Permeability{el2}(b2{sd2});
%                     
%                     % Field values
%                     line1_Bx = 1e3*obj.Flux.x_comp{el1}(b1{sd1});
%                     line1_By = 1e3*obj.Flux.y_comp{el1}(b1{sd1});
%                     
%                     line2_Bx = 1e3*obj.Flux.x_comp{el2}(b2{sd2});
%                     line2_By = 1e3*obj.Flux.y_comp{el2}(b2{sd2});
%                     
% %                     
% %                     line1_B  = obj.Flux.abs{el1}(b1{sd1});
% %                     line2_B  = obj.Flux.abs{el2}(b2{sd2});
% %                     
%                     
%                     % Normal and tangent values
%                     l1_nx = obj.metrics.N{1,el1}(b1{sd1});
%                     l1_ny = obj.metrics.N{2,el1}(b1{sd1});
%                     
%                     l2_nx = obj.metrics.N{1,el2}(b2{sd2});
%                     l2_ny = obj.metrics.N{2,el2}(b2{sd2});
%                     
%                     
%                     % Compute the normal B and tangential H
%                     
%                     l1_Bn = line1_Bx.*l1_nx + line1_By.*l1_ny;
%                     l2_Bn = line2_Bx.*l2_nx + line2_By.*l2_ny;
%                     
%                     l1_Ht = (line1_By.*l1_nx - line1_Bx.*l1_ny)./line1_mu;
%                     l2_Ht = (line2_By.*l2_nx - line2_Bx.*l2_ny)./line2_mu;
%                     
%                      % Computing the average Bn, Ht, B and H
% %                     
%                     l_Bn  = mean([abs(l1_Bn); abs(l2_Bn)]);
%                     
%                     l_Ht  = mean([abs(l1_Ht); abs(l2_Ht)]);
%                     
%                     l1_B    = abs(l_Bn + l_Ht.*line1_mu*1i);
%                     
%                     l2_B    = abs(l_Bn + l_Ht.*line2_mu*1i);
%                                    
%                     obj.Flux.abs{el1}(b1{sd1}) = l1_B;
%                     obj.Flux.abs{el2}(b2{sd2}) = l2_B;
%                 end
%             end
%             
%             % Adjusting the values on the points
%             shared_point = obj.points.shared_point;
%             for k = 1:Np
%                 
%                 % Filter the corner points
%                 if numel(shared_point{1,k}) > 1
%                     % Get the ids of the point k in all shared elements
%                     elid = shared_point{1,k}; lid = shared_point{2,k};
%                     
%                     % Num
%                     els = numel(lid); % Number of shared elements
%                     
%                     % Get the values of the k point on all shared elements
%                     vpx = zeros(1,els); vpy = vpx; vpabs = vpx; vpmu = vpx;
%                     pnx = vpx; pny = vpx;
%                     
%                     pid = zeros(els,4);
%                     
%                     for p = 1:els
%                         eli = elid(p);
%                         
%                         tmp = int_el(obj.Flux.x_comp{eli},'address');
%                         pid(p,:) = tmp{3};
%                         
%                         vpx(p)   = 1e3*obj.Flux.x_comp{eli}(pid(p,lid(p)));
%                         vpy(p)   = 1e3*obj.Flux.y_comp{eli}(pid(p,lid(p)));                         
%                         vpabs(p) =     obj.Flux.abs{eli}(pid(p,lid(p)));
%                         
%                         vpmu(p)  = obj.Permeability{eli}(pid(p,lid(p)));  
%                         
%                         
%                         pnx(p) = obj.metrics.N{1,eli}(pid(p,lid(p)));
%                         pny(p) = obj.metrics.N{2,eli}(pid(p,lid(p)));
%                     end
%                     
%                     
%                     % Compute the normal B and tangential H
%                     
%                     vpbn = vpx.*pnx + vpy.*pny;
%                     vpht = (vpy.*pnx - vpx.*pny)./vpmu;
%                     
%                     vabs = abs(mean(abs(vpbn)) + mean(abs(vpht)).*vpmu*1i);
%                     
%                 end
%             end
                
            
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
                    obj.Potential{k},flux_lines);
            end
        end
        
        %------------------------------------------------------------------
        % Plot a quantity on SEM nodes
        %------------------------------------------------------------------
        function plot_surf_var(obj,var)
            Nel = obj.ProblemData.Nel;
                       
                    for k = 1:Nel
                        surf(obj.mappings.Xm{k},obj.mappings.Ym{k},...
                            var{k},'edgecolor','none');
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
                          lxy        = PolynomialInterp12D(csi_i,eta_i,X,Y);
           
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
        function [Permeability, K ]= ...
                            get_next_permeability(obj,Materials,xmlContent)
            % Initial parameters
            
            options = optimoptions('fsolve','Display','off',...
                'FunctionTolerance',1e-8);
            mod_B = obj.Flux.abs;
            Nel   = obj.ProblemData.Nel;
            mu0   = pi*4e-7;
            Nr    = eval(xml2matlab(xmlContent,'Regions'...
                                                     ,0,'nR','Attribute'));
            parameters = obj.ProblemData.physical_parameters;
            Nparam = length(parameters);
            Permeability = cell(1,Nel);
            K     = cell(1,Nel);
            % Loading in the local workspace
            
            check_list = false(1,Nparam);
            
            for k = 1:Nparam
                try
                    eval([parameters{k,1}, '=',parameters{k,2}, ';']);
                catch
                    check_list(k) = true;
                end
            end
            check_list = find(check_list==true);
            % running the parameters with errors
            
            for k = check_list
                    eval([parameters{k,1}, '=', parameters{k,2}, ';']);
            end
            
            % !!! Temporary analytical BH curve
            bh = @(H) mu0*H+2*Js/pi.*atan((pi*(mur-1)*mu0*H)./(2*Js));
            bhp = @(H) mu0 + (mur-1)*mu0./((pi*(mur-1)*mu0)^2/(2*Js)^2*H.^2+1);
            % Obtaining the magnetic field strength
            mod_H    = cell(1,Nel);
            B_bh = cell(1,Nel);
            
            % Getting the features of each element
            fcn_hdl = false(1,Nel);
            
            for k = 1:Nr
                % Extracting the initial condition
                init_prm  = eval(['[' xml2matlab(obj.xmlContent,...
                    'region',k-1,'InitialPermeability','Attribute') ']']);
                if isempty(init_prm)
                else
                    ElementList  = eval(['[' xml2matlab(obj.xmlContent,...
                    'region',k-1,'ElementList','Attribute') ']']);
                    fcn_hdl(ElementList) = true;
%                     char_fcn = xml2matlab(obj.xmlContent,...
%                                    'region',k-1,'MagneticMaterial','Attribute');
                end
            end
            for k = 1:Nel
                % Computing the magnetic field strength
                mod_H{k} = mod_B{k}./(Materials.Permeability{k}*mu0);
                
                % Computing the remanence of B according to the BH curve
                if fcn_hdl(k)
                    H = mod_H{k};
                    h = fsolve(@(x) bh(x) - mod_B{k}(:),mod_H{k}(:),options);
                    H = reshape(h,size(H));
                    B_bh{k} = bh(H);
                    Brem = B_bh{k}-bhp(H).*H;
                    Permeability{k} = bhp(H); %
                    K{k} = Brem./(mod_B{k}.*bhp(H));
                    
                else
                    B_bh{k} = mod_B{k};
                    Permeability{k} = Materials.Permeability{k};
                    K{k} = zeros(size(mod_B{k}));
                end
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

function [fx, fy] = sem_gradient(f, Lx, Ly, Xxi, Xeta, Yxi, Yeta, J, ...
                                                               ElementSize)
   fxi  = reshape(Lx*f(:),ElementSize(1),ElementSize(2));
                   
   feta = reshape(Ly*f(:),ElementSize(1),ElementSize(2));
                   
                   
   fx = (Yeta.*fxi - Yxi.*feta)./J;
   
   fy = (-Xeta.*fxi + Xxi.*feta)./J;
end