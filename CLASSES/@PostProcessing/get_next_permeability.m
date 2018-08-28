function [Permeability, K, nonlin_elem ]= ...
    get_next_permeability(obj,Materials,xmlContent)
% Initial parameters

%             options = optimoptions('fsolve','Display','off',...
%                 'FunctionTolerance',1e-8);
mod_B = obj.Flux.abs;
% B_x   = obj.Flux.x_comp;
% B_y   = obj.Flux.y_comp;
Nel   = obj.ProblemData.Nel;
%             mu0   = pi*4e-7;
Nr    = eval(xml2matlab(xmlContent,'Regions'...
    ,0,'nR','Attribute'));
parameters = obj.ProblemData.physical_parameters;
Nparam = length(parameters);
Permeability = cell(1,Nel);
% Permeabilitys = cell(1,Nel);
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
%             bh = @(H) mu0*H+2*Js/pi.*atan((pi*(mur-1)*mu0*H)./(2*Js));
%             bhp = @(H) mu0 + (mur-1)*mu0./((pi*(mur-1)*mu0)^2/(2*Js)^2*H.^2+1);
% Obtaining the magnetic field strength
%             mod_H    = cell(1,Nel);
%             B_bh = cell(1,Nel);

% Getting the features of each element
nonlin_elem = false(1,Nel);

for k = 1:Nr
    % Extracting the initial condition
    init_prm  = eval(['[' xml2matlab(obj.xmlContent,...
        'region',k-1,'InitialPermeability','Attribute') ']']);
    if isempty(init_prm)
    else
        ElementList  = eval(['[' xml2matlab(obj.xmlContent,...
            'region',k-1,'ElementList','Attribute') ']']);
        nonlin_elem(ElementList) = true;
        %                     char_fcn = xml2matlab(obj.xmlContent,...
        %                                    'region',k-1,'MagneticMaterial','Attribute');
    end
end

bh = bh_class(1);
for k = 1:Nel
    % Computing the magnetic field strength
    %                 mod_H{k} = mod_B{k}./(Materials.Permeability{k}*mu0);
    
    % Computing the remanence of B according to the BH curve
    if nonlin_elem(k)
        %                     H = mod_H{k};
        %                     h = fsolve(@(x) bh(x) - mod_B{k}(:),mod_H{k}(:),options);
        %                     H = reshape(h,size(H));
        %                     B_bh{k} = bh(H);
        %                     Brem = B_bh{k}-bhp(H).*H;
        %                     Permeability{k} = bhp(H); %
        %                     B_bh{k} = mod_B{k};
        B = mod_B{k};
        if min(B(:)==0) && max(B(:))==0
            B = B+0.01e-3;
        end
        [Brem, ~, Permeability{k}] = bh.BHtool(B);
%     [NU, ~] = BHtoolNR (B*1e3,2);
%         Permeability{k} = 1./NU;
%         Permeabilitys{k} = bh.permeability(B*1e3);
%         Permeability{k} = bh.permeability(B*1e3);
% % % %         Brem = Brem;
        K{k} = Brem./(B.*Permeability{k});
    else
        %                     B_bh{k} = mod_B{k};
        Permeability{k}  = Materials.Permeability{k};
%         Permeabilitys{k} = Materials.Permeability{k};
        K{k} = zeros(size(mod_B{k}));
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