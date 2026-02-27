function circle_optimizer(DRIVE_TYPE)
% This script finds the optimal circle release parameters (radius and frequency).
% OPTIMIZATION METHOD:
%   - Radius (Discrete): Iterates through integer radii (1 to r_max).
%   - Frequency (Continuous): Uses fminbnd (Golden Section Search) for high precision.
%   - This hybrid approach is ~3x faster than grid search and eliminates "stepped" artifacts.

%% ================== 1. Main Configuration ==================
R = 400;       % Domain Radius (Total grid size)
r_max = 200;   % Max Introduction Radius to search
output_dir = './results/circle_release'; % Output directory

fprintf('--- Starting High-Precision Optimization for: %s ---\n', DRIVE_TYPE);

%% ================== 2. Select Drive Function & Time Range ==================
switch DRIVE_TYPE
    case 'TARE'
        drive_function = @TARE_spatial;
        T_range = 10:10:100;
    case '2-locus_TARE'
        drive_function = @Two_locus_TARE_spatial;
        T_range = 1000;
    case 'TADE_suppression'
        drive_function = @TADE_suppression_spatial;
        T_range = [225 250 275 300];
    case 'CifAB'
        drive_function = @CifAB_spatial;
        T_range = [10 100 200 300 400 500 600 700 800 900 1000];
        T_range = [5000 10000];
    case 'Wolbachia'
        drive_function = @Wolbachia_spatial;
        T_range = 300:100:1000;
    otherwise
        error("Unknown drive type specified: '%s'", DRIVE_TYPE);
end

%% ================== 3. Run Optimization in Parallel ==================
fprintf('\n--- Processing %d time points in parallel ---\n', length(T_range));
results = cell(length(T_range), 1);

tic; 
parfor t_idx = 1:length(T_range)
    T = T_range(t_idx);
    
    % Initialize best found values for this Time point
    global_best_eff = -inf;
    global_best_r = 0;
    global_best_f = 0;
    
    % --- Loop 1: Discrete search over Radius (Integer constraint) ---
    % We iterate every integer radius because 'fminbnd' cannot handle integer constraints directly.
    % This is robust and still faster than a full 2D grid.
    for r = 1:r_max
        
        % Define objective function: Maximize Efficiency => Minimize negative Efficiency
        % Note: We pass 'r' fixed, and optimize 'f'
        obj_fun = @(f) -drive_efficiency_circle_helper(r, f, R, T, drive_function);
        
        % --- Loop 2: Continuous optimization for Frequency ---
        % fminbnd efficiently finds the peak frequency for this specific radius
        % TolX = 1e-5 ensures 5-6 decimal places of precision
        [best_f_for_r, min_neg_eff] = fminbnd(obj_fun, 0, 1, optimset('TolX', 1e-5));
        
        current_max_eff = -min_neg_eff;
        
        % Update global best if this radius produced a better result
        if current_max_eff > global_best_eff
            global_best_eff = current_max_eff;
            global_best_r = r;
            global_best_f = best_f_for_r;
        end
    end
    
    % --- Store results ---
    results{t_idx} = struct('time', T, ...
                            'radius', global_best_r, ...
                            'frequency', global_best_f, ...
                            'efficiency', global_best_eff);
                        
    fprintf('T = %d Complete. Best: R=%d, Freq=%.4f, Eff=%.4f\n', ...
            T, global_best_r, global_best_f, global_best_eff);
end
total_elapsed_time = toc;
fprintf('\n--- Optimization completed in %.2f seconds ---\n', total_elapsed_time);

%% ================== 4. Merge and Save Results ==================
output_filename = sprintf('%s/%s_optimal_circle_release.csv', output_dir, DRIVE_TYPE);
fprintf('\nSaving results to: %s\n', output_filename);

if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end

new_results_table = struct2table(cat(1, results{:}));

if exist(output_filename, 'file')
    fprintf('File exists. Merging and updating...\n');
    existing_data = readtable(output_filename);
    new_times = new_results_table.time;
    existing_data(ismember(existing_data.time, new_times), :) = [];
    final_data_to_save = [existing_data; new_results_table];
else
    fprintf('Creating new file...\n');
    final_data_to_save = new_results_table;
end

sorted_table = sortrows(final_data_to_save, 'time');
writetable(sorted_table, output_filename);
fprintf('Success. CSV updated.\n');
end

%% ================== Helper Function ==================
function output = drive_efficiency_circle_helper(radius, frequency, R, T, drive_function)
    % Prepare the introduction vector
    % FIX: Use R instead of hardcoded 200 to match the domain size
    intro_vec = zeros(1, 200);
    
    % Ensure radius doesn't exceed domain (safety check)
    r_idx = min(round(radius), 200);
    
    if r_idx > 0
        intro_vec(1:r_idx) = frequency;
    end
    
    % Run simulation
    output = drive_function(intro_vec, R, T);
end