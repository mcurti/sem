function [ xj, wj ] = LegendreGausLobattoNodesAndWeights1( N )
%LegendreGausLobattoNodesAndWeights Summary of this function goes here
%   Detailed explanation goes here
try
    NodesAndWeights = [];
    load('NodesAndWeights.mat','NodesAndWeights')
    
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

