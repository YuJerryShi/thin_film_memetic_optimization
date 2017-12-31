%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2017 Jerry Shi
% 
% Run this script once to specify the wavelengths of your interest
% This script reads into the dispersion_data folder and interpolates refractive indices 
% The interpolated refractive indices are used to perform the memetic algorithm optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; 

%%
L0 = 1e-6; % Length scale (1e-6 corresponds to microns)
c0 = 2.9979e8 / L0; % Speed of light in vacuum

%% specify the file name
filename = 'ref_ind_cooler'; 
disp_folder = 'dispersion_data'; 
addpath(disp_folder); 

%% Specify the wavelength ranges that you are interested in 
% You can specify different wavelength ranges
% format: [lambda_min1, lambda_max1, Nsamples1;  
%          lambda_min2, lambda_max2, Nsamples2; ...] 
%
% Wavelengths are evenly spaced in frequency

% lambda_vec: coarse wavelength vector used for optimization
% lambda_vec_fine: fine wavelength vector used for plotting

lambda_range = [0.35, 1.80, 100; 
                6.00, 20.0, 150]; 

%% Generate coarse wavelength vector for optimization
N_ranges = length(lambda_range(:, 1)); 
lambda_vec = []; 

for i = 1:N_ranges
    lambda_low_temp = lambda_range(i, 1); 
    lambda_high_temp = lambda_range(i, 2); 
    Nsample_temp = lambda_range(i, 3); 
    
    omega_high_temp = 2*pi*c0 / lambda_low_temp; 
    omega_low_temp = 2*pi*c0 / lambda_high_temp; 
    
    omega_vec_temp = linspace(omega_low_temp, omega_high_temp, Nsample_temp); 
    lambda_vec_temp = 2*pi*c0 ./ omega_vec_temp; 
    
    lambda_vec = [lambda_vec; flipud(transpose(lambda_vec_temp))]; 
end

%% Generate fine wavelength vector for plotting
N_fine = 2000; 
lambda_min = min(lambda_vec); 
lambda_max = max(lambda_vec); 

omega_max = 2*pi*c0 / lambda_min; 
omega_min = 2*pi*c0 / lambda_max; 

omega_vec_fine = linspace(omega_min, omega_max, N_fine); 
lambda_vec_fine = flipud(transpose(2*pi*c0 ./ omega_vec_fine)); 

%%
N = length(lambda_vec); 
num_materials = 10; 

material_key = cell(num_materials, 1); 

ref_ind = zeros(N, num_materials); 
ref_ind_fine = zeros(N_fine, num_materials); 

%% Read and interpolate material refractive indices

% Al2O3 
mat_num = 1; 
material_key{mat_num} = 'Al2O3'; 
A = xlsread('mat_Al2O3'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% HfO2
mat_num = 2; 
material_key{mat_num} = 'HfO2'; 
A = xlsread('mat_HfO2'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% MgF2
mat_num = 3; 
material_key{mat_num} = 'MgF2'; 
A = xlsread('mat_MgF2'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% SiC
mat_num = 4; 
material_key{mat_num} = 'SiC'; 
A = xlsread('mat_SiC'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% SiN
mat_num = 5; 
material_key{mat_num} = 'SiN'; 
A = xlsread('mat_SiN'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% SiO2
mat_num = 6; 
material_key{mat_num} = 'SiO2'; 
A = xlsread('mat_SiO2'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% TiO2
mat_num = 7; 
material_key{mat_num} = 'TiO2'; 
A = xlsread('mat_TiO2'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% Ta2O5
mat_num = 8; 
material_key{mat_num} = 'Ta2O5'; 
A = xlsread('mat_Ta2O5'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% Al
mat_num = 9; 
material_key{mat_num} = 'Al'; 
A = xlsread('mat_Al'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

% Ag
mat_num = 10; 
material_key{mat_num} = 'Ag'; 
A = xlsread('mat_Ag'); 
ref_ind(:, mat_num) = disp_interpolate(lambda_vec, A); 
ref_ind_fine(:, mat_num) = disp_interpolate(lambda_vec_fine, A); 

%%
save(strcat(disp_folder, '/', filename), 'lambda_vec', 'lambda_vec_fine', 'ref_ind', 'ref_ind_fine', 'material_key'); 
rmpath(disp_folder); 