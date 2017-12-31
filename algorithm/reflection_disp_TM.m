function [R_TM] = reflection_disp_TM(n_in, n_out, n_disp, d, lambda, theta_in)
%REFLECTION_DISP_TM Calculates the oblique TE reflection spectrum of a multi-layer structure
%   R_TM: reflection spectrum
% 
%   n_in: scalar/vector that specifies the refractive index of the incident material
%   n_out: scalar/vector that specifies the refractive index of the substrate material
%   n_disp: matrix that specifies the dispersive refractive indices of the structure
%   d: vector that specifies the thicknesses of each layer
%   theta_in: vector that specifies the wavelength of interest

%%
d = abs(d); 
K = length(d); 

%%
theta_out = transpose(asin(n_in * sin(theta_in) / n_out)); 

Z_out = (1./n_out .* cos(theta_out)); 
Z_in = 1/n_in * cos(theta_in); 

%% Iteratively use the impedance method
Z_inter = Z_out; 

for i = 1:K
    j = K-i+1; 
    nj = n_disp(:, j); 
    
    if i == 1
        theta_j = asin(n_out .* sin(theta_out) ./ nj); 
        
    else
        theta_j = asin(n_in .* sin(theta_in) ./ nj); 
    end
    
    dj = d(j); 
    
    Z_inter = 1./nj .* cos(theta_j) .* (Z_inter + 1i*1./nj.*cos(theta_j) .* tan(2*pi*nj./lambda .* cos(theta_j) * dj)) ./ (1./nj.*cos(theta_j) + 1i.*Z_inter .* tan(2*pi*nj./lambda .* cos(theta_j) * dj)); 
    
end

R_TM = abs((Z_inter - Z_in) ./ (Z_inter + Z_in)).^2; 

end

