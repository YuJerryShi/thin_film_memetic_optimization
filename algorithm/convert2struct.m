function [structure] = convert2struct(material_key, n_index_layers, d_layers)
%CONVERT2STRUCT This function populates the ordering of materials and layer
%thicknesses into a readable cell that reprsents the multi-layer structure
% 
%   structure: a 2-column cell that stores the multi-layer structure
%              Left column stores material names, right column stores thicknesses
%   material_key: a cell that contains all the material names
%   n_index_layers: a vector that contains the ordering of materials
%   d_layers: a vector that contains the layer thicknesses

K = length(n_index_layers); 
structure = cell(K, 2); 

for i = 1:K
    structure{i, 1} = material_key{n_index_layers(i)}; 
    structure{i, 2} = d_layers(i); 
end


end

