function [ sx, sy ] = draw_line( varargin )
%draw_line Draws a straight line or a circle arc between two points
%   The function has two modes, sraight line and circle
% In the line mode:
% last argument should be 'line'
% v1 - P1, first point, 2x1 vector for x and y addres
% v2 - P2, second point 2x1 vector for x and y addres
% v3 - arg, csi or eta grid points Nx1 vector
% 
% In the arc mode: !!Note that the center point now is 0 0!!
% The last argument should be 'arc'
% v1 - P1, first point, 2x1 vector for x and y addres
% v2 - P2, second point 2x1 vector for x and y addres
% v3 - arg, csi or eta grid points Nx1 vector
mode = varargin{nargin};

if strcmp(mode,'segment')
   mode = 'line';
end

n2range = @(x,a1,a2) (x-min(x))./(max(x)-min(x))*(a2-a1)+a1;
x = 1; y = 2;
switch mode
    case 'line'
        P1  = varargin{1};
        P2  = varargin{2};
        arg = varargin{3};
        sx  = n2range(arg,P1(x),P2(x))';
        sy  = n2range(arg,P1(y),P2(y))';
%         if (P2(y)-P1(y))==0
%             sx  = x1;
%             sy  = ones(1,length(y1))*P1(y);
%         elseif (P2(x)-P1(x))==0
%             sx  = ones(1,length(y1))*P1(x);
%             sy  = y1;
%         else
%             sx  = (y1 - P1(y))/(P2(y)-P1(y))*(P2(x)-P1(x))+P1(x);
%             sy  = (x1 - P1(x))/(P2(x)-P1(x))*(P2(y)-P1(y))+P1(y);
%         end
    case 'arc'
        P1  = varargin{1};
        P2  = varargin{2};
        arg = varargin{3};
        a2     = n2range(arg,cart2pol(P1(x),P1(y)),cart2pol(P2(x),P2(y)))';
        [~,R]  = cart2pol(P1(x),P1(y)); [~,Rt]  = cart2pol(P2(x),P2(y));
        if abs(round(R-Rt,10))>0
            error('Points are not concentric!!!')
        end
        sx     = R*cos(a2);
        sy     = R*sin(a2);
    case 'arcc'
        P1  = varargin{1};
        P2  = varargin{2};
        arg = varargin{3};
        Po = varargin{4};
        %------------------------------------------------------------------
        R = abs(abs(P1(x) - Po(x))+abs(P1(y) - Po(y))*1i);
        alpha_1 = atan2(P1(y)-Po(y),P1(x)-Po(x));
        alpha_2 = atan2(P2(y)-Po(y),P2(x)-Po(x));
%         if alpha_2< -pi/2
%             alpha_2 = 2*pi - abs(alpha_2);
% %             disp('1')
%         end
        if alpha_1 < 0
            alpha_1 = 2*pi + alpha_1;
        end
        
        if alpha_2 < 0
            alpha_2 = 2*pi + alpha_2;
        end
        %------------------------------------------------------------------
        a2     = n2range(arg,alpha_1,alpha_2);
        sx     = Po(x) + R*cos(a2)';
        sy     = Po(y) + R*sin(a2)';
    case 'arcangle'
        P1    = varargin{1};
        P2    = varargin{2};
        arg   = varargin{3};
        alpha = varargin{4};
        cx  = (P1(x)+P2(x))*.5; cy  = (P1(y)+P2(y))*.5;
        snx   = sign(P1(x) - cx);
        sny   = sign(P1(y) - cy);
        slope   = -(P2(x)-P1(x))/(P2(y)-P1(y));
        alpha_1 = abs(atan(slope));
        beta    = alpha_1 - .5*alpha;
        
        
        B2 = abs(abs(P1(x) - P2(x))+abs(P1(y) - P2(y))*1i);
        R  = B2/(2*sin(alpha/2));
        
        offx = abs(R*cos(beta));
        offy = abs(R*sin(beta));
        if snx == -1 && sny == 0
        Po(1) = P1(1) + offx;
        Po(2) = P1(2) - offy;
        elseif snx == 0 && sny == 1
        Po(1) = P1(1) - offx;
        Po(2) = P1(2) - offy;
        elseif snx == 1 && sny == 0
        Po(1) = P1(1) - offx;
        Po(2) = P1(2) + offy;
        elseif snx == 0 && sny == -1
        Po(1) = P1(1) + offx;
        Po(2) = P1(2) + offy;
        else
        Po(1) = P1(1) - offx;
        Po(2) = P1(2) - offy;
        end
%         Po
        angle_1 = cart2pol(P1(x)-Po(x),P1(y)-Po(y));
        angle_2 = cart2pol(P2(x)-Po(x),P2(y)-Po(y));
        if angle_1< -pi/2
            angle_1 = 2*pi - abs(angle_1);
%             disp('1')
        end
%         if angle_2< -pi/2
%             angle_2 = 2*pi - abs(angle_2);
%             disp('2')
%         end
        a2      = n2range(arg,angle_1,angle_2);
        
        sx     = Po(x) + R*cos(a2);
        sy     = Po(y) + R*sin(a2);
    otherwise
        error('Input mode not recognized')
end

end

