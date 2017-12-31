function [n_material_layers, d_layers] = generate_population(num_population, num_layers, list_materials, d_max)
%GENERATE_POPULATION Randomly generates the population used in the memetic
%algorithm
% 
%   n_material_layers: Matrix that stores the ordering of the materials in the population
%                      size: [num_population, num_layers]
%   d_layers: Matrix that stores the ordering of layer thicknesses in the population 
%             size: [num_population, num_layers]
% 
%   num_population: Integer that specifies the number of population
%   num_layers: Integer that speficies the number of layers
%   list_materials: Vector of integers that specifies which materials are used 
%   d_max: Continuous variable that limits the maximum layer thickness

max_index = length(list_materials); 
temp_ordering = randi(max_index, num_population, num_layers); 


%% Make sure there's no duplicates
for i = 1 : num_population
    for j = 2 : num_layers-1
        if temp_ordering(i, j) == temp_ordering(i, j-1) && temp_ordering(i, j) == temp_ordering(i, j+1)
            Z = 1:max_index; 
            Z(Z == temp_ordering(i, j-1)) = []; 
            Z(Z == temp_ordering(i, j+1)) = []; 
            
            temp_ordering(i, j) = datasample(Z, 1); 
            
        elseif temp_ordering(i, j) == temp_ordering(i, j-1)
            Z = 1:max_index; 
            Z(Z == temp_ordering(i, j-1)) = []; 
            
            temp_ordering(i, j) = datasample(Z, 1); 
            
        elseif temp_ordering(i, j) == temp_ordering(i, j+1)
            Z = 1:max_index; 
            Z(Z == temp_ordering(i, j+1)) = []; 
            
            temp_ordering(i, j) = datasample(Z, 1); 
            
        end
    end
end

%% Convert temp_index to n_layers_index
n_material_layers = list_materials(temp_ordering); 

%% Genearte layer thicknesses
d_layers = d_max/2*ones(num_population, num_layers) + d_max * (rand(num_population, num_layers)-0.5); 

end



