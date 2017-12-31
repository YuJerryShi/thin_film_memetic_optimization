%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2017 Jerry Shi
% 
% Run this script to perform the memetic algorithm optimization
% You should run disp_generate.m once first to specify the wavelengths of your interest. 
% Below, you can specify your objective function
% You can also tune the parameters of the memetic algorithm
% The optimized structure will be shown in the Command Window
% 
% Details of the algorithm: see article DOI 10.1021/acsphotonics.7b01136
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; 

%% Optional: Seed the random number generator for reproducible optimization runs
rng(1, 'twister') % Comment this out if you want different optimization results each time

%% Graphics parameters
lw = 2; % Linewidth of curves
ft = 12; % Font size
addpath('algorithm'); 

%% Load the refractive index data
L0 = 1e-6; % Length scale (1e-6 corresponds to microns)
load dispersion_data/ref_ind_cooler

%% Wavelengths
N = length(lambda_vec); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the target reflection spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target
R_target = 0*((lambda_vec > 8) & (lambda_vec < 13)) + ...
         + 1*(lambda_vec < 8) + ...
         + 1*(lambda_vec > 13);  
     
% Weight of each wavelength
W = eye(N); 
for i = 1:N
    if lambda_vec(i) < 4
        W(i, i) = 5; 
    elseif lambda_vec(i) > 13
        W(i, i) = 1; 
    else
        W(i, i) = 1; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the materials used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material indices: 
% 1. Al2O3
% 2. HfO2
% 3. MgF2
% 4. SiC
% 5. SiN
% 6. SiO2
% 7. TiO2
% 8. Ta2O5
% 9. Al
% 10. Ag

list_materials = [1 2 5 6 7 8]; % Use materials Al2O3, HfO2, SiN, SiO2, TiO2, Ta2O5

% Incidence material
n_in = 1; % Incident from air. 

% Substrate material
n_material_sub = 10; % Ag substrate
n_sub = ref_ind(:, n_material_sub);
n_sub_fine = ref_ind_fine(:, n_material_sub); 


% % Optional % % 
% Fixed materials between substrate and layered structure
% Use n_material_fixed_sub = []; d_fixed_sub = []; if you don't want this
n_material_fixed_sub = []; 
d_fixed_sub = []; 

%% Definte the angle of incidence. Note: nonzero theta is usually slow
theta_in = 0; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Memetic algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for generating population and maximum number of generations
M = 1200; % Population size
K = 6; % Number of layers
d_max = 0.80; % Maximum thickness of each layer

max_iter = 61; % Maximum number of generations

%% Mutation rate
mutation_rate = 0.05; % A typical value is between 0.01 and 0.15

%% Reselection rates
% Rule: 
% Keep top RATE_TOP of the population as parents 
% Select the rest of the RATE_GOOD of the parents from the top-half of population
% Select the remaining RATE_POOR of the parents from the bottom half

rate_top = 0.00; % A typical value is between 0 and 0.02
rate_good = 0.85; % A typical value is between 0.80 and 0.99
rate_poor = 1 - rate_good; 

%% Refinement rate
refine_num = 60; % Number of individuals that are refined
refine_period = 20; % Refinement is done for every "refine_period" generations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the memetic algorithm main program located in "/algorithm"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
memetic_algorithm_main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the best structure's spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_best = ref_ind(:, n_material_best); 
n_best_fine = ref_ind_fine(:, n_material_best); 

% Display the best structure (including the fixed materials above substrate)
structure_best = convert2struct(material_key, [n_material_best, n_material_fixed_sub], [d_best, d_fixed_sub])

% Compute the final optimized spectrum and the spectrum without multi-layers
R_best_fine = reflection_disp(n_in, n_sub_fine, [n_best_fine, n_fixed_sub_fine], [d_best, d_fixed_sub], lambda_vec_fine, theta_in); 
R_no_layers_fine = reflection_disp(n_in, n_sub_fine, n_fixed_sub_fine, d_fixed_sub, lambda_vec_fine, theta_in); 

% Plotting the optimized spectrum
figure; hold on; 
scatter(lambda_vec, 1-R_target, 'go'); 
plot(lambda_vec_fine, 1-R_best_fine, 'b', 'linewidth', lw); 
plot(lambda_vec_fine, 1-R_no_layers_fine, 'r--', 'linewidth', lw); 
hold off; 
axis([0.3 20 0 1])
xlabel('wavelength (\mum)'); ylabel('Emissivity'); 
legend('Target', 'Optimum structure', 'Without multi-layer'); 
set(gca, 'fontsize', ft); 
box on; 


