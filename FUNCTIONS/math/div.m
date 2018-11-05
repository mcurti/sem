function out = div(u,v)
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
        % Option for vector vector
        % [u1 u2]*[v1 v2]' = u1*v1 + u2*v2
        out = u(:,:,1)*v(:,:,1) + u(:,:,2)*v(:,:,2);
    case 'tensorvector'
        % Option for tensor vector
        % [u1 u3; u2 u4]*[v1 v2]' = [u1*v1 + u3*v2; u2*v1 + u4*v2]
        out = zeros(size(v));
        out(:,:,1) = u(:,:,1)*v(:,:,1) + u(:,:,3)*v(:,:,2);
        out(:,:,2) = u(:,:,2)*v(:,:,1) + u(:,:,4)*v(:,:,2);
        
end
end

