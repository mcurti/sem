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

check_list = false(1,K);

for k = 1:K
    try
        eval([parameters{k,1}, '=',parameters{k,2}, ';']);
        assignin('caller',parameters{k,1},eval(parameters{k,1}))
    catch
        check_list(k) = true;
    end
end
check_list = find(check_list==true);
% running the parameters with errors
check_list1 = false(1,K);
for k = check_list
    try
        eval([parameters{k,1}, '=',parameters{k,2}, ';']);
        assignin('caller',parameters{k,1},eval(parameters{k,1}))
    catch
        check_list1(k) = true;
    end
end

check_list1 = find(check_list1==true);
for k = check_list1
    eval([parameters{k,1}, '=',parameters{k,2}, ';']);
    assignin('caller',parameters{k,1},eval(parameters{k,1}))
end

end

