function [ref_ind] = disp_interpolate(lambda_in, A)
%DISP_INTERPOLATE Summary of this function goes here
%   Detailed explanation goes here

N = length(lambda_in); 
N2 = length(A(:, 1)) ; 

ref_ind = zeros(N, 1); 

lambda0 = A(:, 1); 
n0 = A(:, 2); 
k0 = A(:, 3); 

%%

for i = 1:N
    for j = 1:N2
        if lambda_in(i) > lambda0(j) && lambda_in(i) <= lambda0(j+1)
%             ref_ind(i) = A(j, 2) - 1i*A(j, 3); 
            ref_ind(i) = n0(j) + (n0(j+1) - n0(j))/(lambda0(j+1) - lambda0(j)) * (lambda_in(i)-lambda0(j)) + ...
                       -1i * (k0(j) + (k0(j+1) - k0(j))/(lambda0(j+1) - lambda0(j)) * (lambda_in(i)-lambda0(j))); 
            break; 
        end
    end
end

end

