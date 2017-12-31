function [n_material_after, d_after] = mutation(n_material_before, d_before, list_materials, d_max)
%MUTATION Randomly perturbs the property of a random layer
%   n_material_after: ordering of the materials after mutation
%   d_after: thicknesses of the materials after mutation
% 
%   n_material_before: ordering of the materials before mutation
%   d_before: thicknesses of the materials before mutation
%   list_materials: list of possible materials to mutate into 
%   d_max: maximum thickness of the mutated layer

%%
K = length(d_before); 
d_after = d_before; 
n_material_after = n_material_before; 

%%
layer_change = randi(K); 

d_after(layer_change) = d_max/2 + d_max * (rand()-0.5); 

index_max = length(list_materials); 
temp_index = randi(index_max); 
n_material_after(layer_change) = list_materials(temp_index); 

end

