function [R] = reflection_disp(n_in, n_out, n_disp, d, lambda, theta_in)
%REFLECTION_DISP Calculates the reflection spectrum of a multi-layer structure
%   R: reflection spectrum
% 
%   n_in: scalar/vector that specifies the refractive index of the incident material
%   n_out: scalar/vector that specifies the refractive index of the substrate material
%   n_disp: matrix that specifies the dispersive refractive indices of the structure
%   d: vector that specifies the thicknesses of each layer
%   theta_in: vector that specifies the wavelength of interest

if theta_in == 0
    %% Normal reflection
    R = reflection_disp_norm(n_in, n_out, n_disp, d, lambda); 
else
    %% Oblique reflection
    R_TE = reflection_disp_TE(n_in, n_out, n_disp, d, lambda, theta_in); 
    R_TM = reflection_disp_TM(n_in, n_out, n_disp, d, lambda, theta_in); 
    
    R = (R_TE + R_TM) / 2; 
end

end

