function out = div(u,v,tu,tv)
%GRAD executes the divergence between  two twodimmensional vectors each
%component of these vectors are expected to be matrices which have to have
%the same size.


sizeu = size(u); sizev = size(v);

if sizeu(3) == 2 && sizev(3) == 2
    option = 'vector';
elseif sizeu(3) == 4 && sizev(3) == 2
    option = 'tensorvector';
end

switch option
    case 'vector'
        %------------------------------------------------------------------
        u1 = u(:,:,1); u2 = u(:,:,2);
        v1 = v(:,:,1); v2 = v(:,:,2);
        %------------------------------------------------------------------
        % Option for vector vector
        % [u1 u2]*[v1 v2]' = u1*v1 + u2*v2
        if tu == 1 && tv == 1
            out = u1*v1 + u2*v2;
        else
            out = mprod(u1,tu,v1,tv) + mprod(u2,tu,v2,tv);
        end
    case 'tensorvector'
        %------------------------------------------------------------------
        u1 = u(:,:,1); u2 = u(:,:,2); u3 = u(:,:,3); u4 = u(:,:,4);
        v1 = v(:,:,1); v2 = v(:,:,2); % v3 = u(:,:,3); v4 = u(:,:,4);
        %------------------------------------------------------------------
        % Option for tensor vector
        % [u1 u3; u2 u4]*[v1 v2]' = [u1*v1 + u3*v2; u2*v1 + u4*v2]
        out = zeros(size(v));
        if tu == 1 && tv == 1
            out(:,:,1) = u1*v1 + u3*v2;
            out(:,:,2) = u2*v1 + u4*v2;
        else
            out(:,:,1) = mprod(u1,tu,v1,tv) + mprod(u3,tu,v2,tv);
            out(:,:,2) = mprod(u2,tu,v1,tv) + mprod(u4,tu,v2,tv);
        end
        
end
end

