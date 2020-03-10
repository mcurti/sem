function figure_config( varargin )
% Function which will make the same configuration for all figure.
% The order of input arguments:
% 1. The figure number
% 2. The figure width
% 3. The figure height


% Definition of the default settings
f_number = 1; w = 10; h = 10; fontsize = 6;

% Altering the settings according to inputs
switch nargin
    case 1
        f_number  = varargin{1}; 
    case 2
        f_number  = varargin{1}; w = varargin{2};
    case 3
        f_number  = varargin{1}; w = varargin{2}; h = varargin{3};
    case 4
        f_number  = varargin{1}; w = varargin{2}; h = varargin{3};
        fontsize = varargin{4};
end
         

f = figure(f_number);

% Units of the figure
f.PaperUnits = 'centimeters';
f.Units      = 'centimeters';
% Capturing the position
pos = get(gcf,'Position');
pos_x = pos(1); pos_y = pos(2);
f.PaperSize = [w h];
f.PaperPosition = [pos_x pos_y w h];
f.Position = [pos_x pos_y w h];
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',fontsize);
end

