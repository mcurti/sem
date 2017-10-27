function [ xj, wj ] = LegendreGausLobattoNodesAndWeights( N )
%LegendreGausLobattoNodesAndWeights Summary of this function goes here
%   Detailed explanation goes here
try
    load('NodesAndWeights.mat')
    
    xj = NodesAndWeights.nodes{N};
    wj = NodesAndWeights.weights{N};
    
    if isempty(xj)
        error('No input for this N')
    end
catch
    syms x
    y  = legendreP(N+1,x) - legendreP(N-1,x);
    
    roots = vpasolve(y == 0);
    
    xj = double(roots);
    
    wj = 2./(N*(N+1)*legendreP(N,xj).^2);
    
    NodesAndWeights.nodes{N}   = xj;
    NodesAndWeights.weights{N} = wj;
    
    save('NodesAndWeights.mat','NodesAndWeights');
end
end

