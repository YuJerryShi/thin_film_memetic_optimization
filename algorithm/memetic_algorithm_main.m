%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs the memetic algorithm given pre-specified parameters
% It displays the progression of the Merit Function
% It outputs the best structure at a given generation in the Command Window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate population (Size M, K layers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n_material_layers, d_layers] = generate_population(M, K, list_materials, d_max); % 

%% Define the merit function 
% Take into account of the fixed materials above the substrate
n_fixed_sub = ref_ind(:, n_material_fixed_sub); 
n_fixed_sub_fine = ref_ind_fine(:, n_material_fixed_sub); 

MF = @(n_layers, d) norm(W * (reflection_disp(n_in, n_sub, [n_layers, n_fixed_sub], abs([d, d_fixed_sub]), lambda_vec, theta_in) - R_target)).^2; 

% Reference fitness function without multi-layers
fitness_no_layers = MF([], []); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize arrays to keep track of population development
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fitness_best = zeros(max_iter, 1); % Keeps track of the fittest individual for each generation
fitness_avg = zeros(max_iter, 1); % Keeps track of the population's average fitness for each generation
fitness_best_so_far = 1e9; % Keeps track of the best fitness so far

fitness_population = zeros(2*M, max_iter); 
fitness_parents = zeros(M, 1); 

d_parents = d_layers; 
n_material_parents = n_material_layers; 

d_childern = zeros(M, K); 
n_material_childern = zeros(M, K); 

%% Calculate how the parents are chosen for the next generation
% Keep top RATE_TOP of the population as parents 
% Select the rest of the RATE_GOOD of the parents from the top-half of population
% Select the remaining RATE_POOR of the parents from the bottom half

N_top = round(M * rate_top); 
N_good = round(M * (rate_good - rate_top)); 
N_poor = round(M * (rate_poor)); 

N_good_pool = M - N_top; 
N_poor_pool = M; 

%% Convergence curve
figure(1); 
set(gca, 'xlim', [0 max_iter+1], 'fontsize', ft); 
xlabel('Generation'); ylabel('log_1_0((Merit Function))'); 
line([0 max_iter+1], log10([fitness_no_layers fitness_no_layers]), 'linestyle', '--', 'color', 'g', 'linewidth', lw); 
box on; 

%% Memetic algorithm iterations
for i = 1:max_iter
    
    %% One point crossover 
    % Shuffle the layers
    shuffle_index = randperm(M); 
    d_parents = d_parents(shuffle_index, :); 
    n_material_parents = n_material_parents(shuffle_index, :); 
    
    % Select crossover point and crossover
    for j = 1:round(M/2)
        d_parent1 = d_parents(2*j-1, :); 
        d_parent2 = d_parents(2*j, :); 
        n_material_parent1 = n_material_parents(2*j-1, :); 
        n_material_parent2 = n_material_parents(2*j, :); 
        
        
        % Choose crossover point
        cross_pt = randi(K-1); 
        
        d_childern(2*j-1, :) = [d_parent1(1:cross_pt) d_parent2(cross_pt+1:end)]; 
        n_material_childern(2*j-1, :) = [n_material_parent1(1:cross_pt) n_material_parent2(cross_pt+1:end)]; 
        
        d_childern(2*j, :) = [d_parent2(1:cross_pt) d_parent1(cross_pt+1:end)]; 
        n_material_childern(2*j, :) = [n_material_parent2(1:cross_pt) n_material_parent1(cross_pt+1:end)]; 
    end
    
    %% Mutate with a moderate probability
    for m = 1:M
        % Mutate all childern with a probability of mutation_rate
        if rand() < mutation_rate
            [n_material_childern(m, :), d_childern(m, :)] = mutation(n_material_childern(m, :), d_childern(m, :), list_materials, d_max); 
        end
    end

    %% Evaluate fitness
    d_generation = [d_parents; d_childern]; 
    n_material_generation = [n_material_parents; n_material_childern]; 
    
    for m = 1:2*M
        n_generation = ref_ind(:, n_material_generation(m, :)); 
        fitness_population(m, i) = MF(n_generation, d_generation(m, :)); 
    end
    
    % Sort fitnesses
    [fitness_arr, fitness_ind] = sort(fitness_population(:, i));   
    
    
    %% Record the history of the fitness function
    fitness_best(i) = fitness_arr(1); 
    fitness_avg(i) = mean(fitness_arr(1:M)); 
    
    % Plot the convergence curve in real time
    figure(1)
    hold on;
    if (i >= 2)
        line([i-1 i], [log10(fitness_best(i-1)) log10(fitness_best(i))], 'color', 'b', 'linewidth', lw); 
        line([i-1 i], [log10(fitness_avg(i-1)) log10(fitness_avg(i))], 'color', 'r', 'linewidth', lw); 
    end
    
    if (i == 2)
        legend('Merit function with no multi-layer', 'Best merit function', 'Average merit function'); 
    end
    hold off; 
    pause(0.01); 
    
    
    %% Hold onto the best device so far
    if fitness_best(i) < fitness_best_so_far
        i
        n_material_best = n_material_generation(fitness_ind(1), :); 
        d_best = d_generation(fitness_ind(1), :);
        
        % Displays the best structure so far and its fitness
        fitness_best_so_far = fitness_best(i)
        structure_best_so_far = convert2struct(material_key, n_material_best, d_best)
    end
    
    %% Update parents arrays according to fitness
    % Keep the N_top best structures
    selected_ind_top = fitness_ind(1 : N_top); 
        
    % Good structures
    pool_good = fitness_ind(N_top+1 : N_top + N_good_pool); 
    temp_rand_index = randperm(N_good_pool) + N_top; 
    selected_ind_good = fitness_ind(temp_rand_index(1 : N_good)); 
        
    % Poor structures
    pool_bad = fitness_ind(N_top + N_good_pool+1 : end); 
    
    temp_rand_index = randperm(N_poor_pool) + N_top+N_good; 
    
    selected_ind_poor = fitness_ind(temp_rand_index(1 : N_poor)); 
    
    % Concatenate everything
    selected_ind = [selected_ind_top; selected_ind_good; selected_ind_poor]; 
    
    d_parents = d_generation(selected_ind, :); 
    
    n_material_parents = n_material_generation(selected_ind, :); 
    
    
        
    %% Refine elites for every m generations
    if mod(i, refine_period) == 0
%     if 1
        %% Re-sort fitness
        for m = 1:M
            n_parents = ref_ind(:, n_material_parents(m, :)); 
            fitness_parents(m) = MF(n_parents, d_parents(m, :)); 
        end

        [fitness_parents, fitness_parents_ind] = sort(fitness_parents); 

        n_material_parents = n_material_parents(fitness_parents_ind, :); 
        d_parents = d_parents(fitness_parents_ind, :); 
        
        %% Local optimization
        for m = 1:refine_num
            m; 
            n_layers_m = ref_ind(:, n_material_parents(m, :)); 

            obj_fun = @(x) MF(n_layers_m, x); 

            x0 = d_parents(m, :); 

            options = optimoptions('fminunc'); 
            options.TolFun = 1e-8; 
            options.TolX = 1e-8; 
            options.MaxFunEvals = 1e5; 
            options.MaxIter = 1e5; 
            options.Display = 'off'; 
            options.Algorithm = 'quasi-newton'; 
            
            [x, ~, ~, output] = fminunc(obj_fun, x0, options); 
            output; 
            
            d_parents(m, :) = x; 
            obj_fun(x); 
        end
    end
    
end

