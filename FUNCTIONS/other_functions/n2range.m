function [ nx ] = n2range( x, a1, a2 )
%n2range scalling a vector from its min and max values to the input a1 and
% a2 numbers
[rs, cs] = size(x);

nx = (x(:)-min(x(:)))./(max(x(:))-min(x(:)))*(a2-a1)+a1;

nx = reshape(nx,rs,cs);
end

