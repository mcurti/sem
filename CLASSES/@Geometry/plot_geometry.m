function plot_geometry(obj,varargin)
mode = varargin{nargin-1};
switch mode
    % Plot points
    case 'points'
        PointsNumber = eval(xml2matlab(obj.GeometryElement...
            ,'Points',0,'p','Attribute'));
        if nargin > 2
            PropertiesNumber = (nargin-2)/2;
        end
        for k = 1:PointsNumber
            hp = plot(obj.points.coordinates(k,1),...
                obj.points.coordinates(k,2),'.');
            %                          text(obj.points.coordinates(k,1)-1,obj.points.coordinates(k,2),sprintf('%d',k),'HorizontalAlignment','right','Color','red')
            
            if exist('PropertiesNumber','var')
                for ii = 1:PropertiesNumber
                    set(hp,varargin{ii*2-1},varargin{ii*2})
                end
            end
        end
        % Plot lines
        
    case 'lines'
        LinesNumber = eval(xml2matlab(obj.GeometryElement...
            ,'Lines',0,'l','Attribute'));
        if nargin > 2
            PropertiesNumber = (nargin-2)/2;
        end
        for k = 1:LinesNumber
            temp_line = obj.lines.vector{k};
            hp = plot(temp_line(1,:),temp_line(2,:));
                               middle = round(numel(temp_line(1,:))/2);
%                                text(temp_line(1,middle)-1,temp_line(2,middle)-0.05,sprintf('%d',k),'HorizontalAlignment','right','Color','green')
            
            if exist('PropertiesNumber','var')
                for ii = 1:PropertiesNumber
                    set(hp,varargin{ii*2-1},varargin{ii*2})
                end
            end
        end
        % Plot element grid
        
    case 'ElementGrid'
        ElementNumber = ...
            eval(xml2matlab(obj.GeometryElement,...
            'Elements',0,'el','Attribute'));
        if nargin > 2
            PropertiesNumber = (nargin-2)/2;
        end
        for k = 1:ElementNumber
            hp = plot(obj.mappings.Xm{k},obj.mappings.Ym{k},'.');
            
            if exist('PropertiesNumber','var')
                for ii = 1:PropertiesNumber
                    set(hp,varargin{ii*2-1},varargin{ii*2})
                end
            end
        end
        
    case 'Element_id'
        ElementNumber = ...
            eval(xml2matlab(obj.GeometryElement,...
            'Elements',0,'el','Attribute'));
        
        for k = 1:ElementNumber
            
            meanx = mean(obj.mappings.Xm{k}(:));
            meany = mean(obj.mappings.Ym{k}(:));
            text(meanx,meany,sprintf('%d',k),'Color','black')
            
            
        end
    case 'Line_id'
        LinesNumber = eval(xml2matlab(obj.GeometryElement...
            ,'Lines',0,'l','Attribute'));
        
        for k = 1:LinesNumber
            
            middle = mean(obj.lines.vector{k},2);
            text(middle(1),middle(2),sprintf('%d',k),'HorizontalAlignment','right','Color','blue')
            
            
        end
    case 'Point_id'
        PointsNumber = eval(xml2matlab(obj.GeometryElement...
            ,'Points',0,'p','Attribute'));
        
        for k = 1:PointsNumber
            
            text(obj.points.coordinates(k,1),obj.points.coordinates(k,2),sprintf('%d',k),'HorizontalAlignment','left','Color','red')
            
            
        end
        
end
end