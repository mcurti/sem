function [FL, remFL] = flux_linkage(obj,elements)

FL = 0; remFL = 0;
for k = elements
    
    S = sum(obj.metrics.W{k}(:).*obj.metrics.J{k}(:));
    
    FL = FL + sum(obj.Potential{k}(:).*obj.metrics.W{k}(:).*obj.metrics.J{k}(:))/S;
    
    
    
    remFL = remFL + sum(obj.remPotential{k}(:).*...
        obj.metrics.W{k}(:).*obj.metrics.J{k}(:))/S;
end

end