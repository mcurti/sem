%------------------------------------------------------------------
% Function to interpolate on FEM nodes, used to check the
% correctness of the problem
%------------------------------------------------------------------
function obj = interpolate_on_FEM(obj,filenames,no_elements)

FluxData  = Flux2DmeshData(filenames);
FEM       = FluxData.FEMfilesInOne;
try
    FEM.noelements = no_elements;
catch
end
sem = 0;
fileName = sprintf('lxy%d.mat',numel(FEM.all_x_nodes));
if exist(fileName,'file')
    load(fileName)
    FEM.int  = sem.int;
    obj.SEM_interpolation.solutions = ...
                                obj.PProcessing.InterpolateOnFEMNodes(FEM);
    obj.SEM_interpolation.FEM       = FEM.all_values;
    obj.SEM_interpolation.TRI = FEM.TRI;
    obj.SEM_interpolation.all_x_nodes = FEM.all_x_nodes;
    obj.SEM_interpolation.all_y_nodes = FEM.all_y_nodes;
else
    sem = obj.PProcessing.InterpolateOnFEMNodes(FEM);
    
    obj.SEM_interpolation.solutions = sem;
    obj.SEM_interpolation.FEM       = FEM.all_values;
    obj.SEM_interpolation.TRI = FEM.TRI;
    obj.SEM_interpolation.all_x_nodes = FEM.all_x_nodes;
    obj.SEM_interpolation.all_y_nodes = FEM.all_y_nodes;
    save(fileName,'sem');
end

end