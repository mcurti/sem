%==========================================================================
% parameter_loader.m
% Created: 09.02.2018 - 09:59:23
% By: M. Curti
%
% Loads the parameters from a structure which contains the variable name
% and value extracted from .xml file
%==========================================================================
function parameter_loader( parameters )
[K, ~] = size(parameters);
% Loading in the local workspace

check_list = true(1,K);
prev_list  = true(1,K);
while sum(check_list)>0
    for k = 1:K
        try
            eval([parameters{k,1}, '=',parameters{k,2}, ';']);
            assignin('caller',parameters{k,1},eval(parameters{k,1}))
            check_list(k) = false;
        catch
            check_list(k) = true;
        end
    end
    difference = xor(check_list,prev_list);
    
    if sum(difference)==0
        
        error(...
        'An error occured while definig the parameters %s\n'...
        ,parameters{check_list,1});
    else
        prev_list = check_list;
    end
end

end

