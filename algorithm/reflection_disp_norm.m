function [R] = reflection_disp_norm(n_in, n_out, n_disp, d, lambda)
%REFLECTION_DISP_NORM Calculates the normal reflection spectrum of a multi-layer structure
%   R: reflection spectrum
% 
%   n_in: scalar/vector that specifies the refractive index of the incident material
%   n_out: scalar/vector that specifies the refractive index of the substrate material
%   n_disp: matrix that specifies the dispersive refractive indices of the structure
%   d: vector that specifies the thicknesses of each layer
%   theta_in: vector that specifies the wavelength of interest

%%
d = abs(d); 
K = length(d); 

Z_out = 1./n_out; 
Z_in = 1/n_in; 

%% Iteratively use the impedance method
Z_inter = Z_out; 

for i = 1:K
    j = K-i+1; 
    nj = n_disp(:, j); 
    dj = d(j); 
    
    Z_inter = 1./nj .* (Z_inter + 1i*1./nj .* tan(2*pi*nj./lambda * dj)) ./ (1./nj + 1i.*Z_inter .* tan(2*pi*nj./lambda * dj)); 
    
end

R = abs((Z_inter - Z_in) ./ (Z_inter + Z_in)).^2; 

end

