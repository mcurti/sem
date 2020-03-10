%==========================================================================
% Physics.m
% Created: 27.11.2017 - 09:07:44
% By: M. Curti
%
% Change a parameter in the workspace of a class
%==========================================================================
function obj = change_parameter( obj,Parameter,Value )

% Estimating the number of parameters
Nparam = length(obj.parameters);

% For loop to change the input parameter
for k = 1:Nparam
    if strcmp(obj.parameters{k,1},Parameter)
        obj.parameters{k,2} = Value;
        return
    end
end
error('Requested parameter does not exist')
end

