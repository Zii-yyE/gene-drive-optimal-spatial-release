function run_circle_optimizer_tade_suppression()
    % Wrapper script for command-line execution for CifAB
    try
        % ================== 1. Main Configuration ==================
        DRIVE_TYPE = 'TADE_suppression';
        R = 400;       % Domain Radius (Total grid size)
        r_max = 200;   % Max Introduction Radius to search
        output_dir = './results/circle_release'; % Output directory
        
        fprintf('--- Starting High-Precision Optimization for: %s ---\n', DRIVE_TYPE);
        
        % ================== 2. Setup ==================
        % Explicitly define the drive function
        drive_function = @TADE_suppression_spatial;
        
        % Time Range (Using the larger range from your snippet)
        % If you want the smaller range, uncomment the line below:
        T_range = [5000 10000 20000];
        
        % ================== 3. Parallel Pool Setup ==================
        % Cleanup old pools to avoid hanging
        delete(gcp('nocreate')); 
        parpool(); % Start new pool using available cores
        
        fprintf('\n--- Processing %d time points in parallel ---\n', length(T_range));
        results = cell(length(T_range), 1);
        
        tic;
        parfor t_idx = 1:length(T_range)
            % [CRITICAL] Prevent worker thread contention
            maxNumCompThreads(1);
            
            T = T_range(t_idx);
            
            % Initialize best found values for this Time point
            global_best_eff = -inf;
            global_best_r = 0;
            global_best_f = 0;
            
            % --- Loop 1: Discrete search over Radius ---
            for r = 1:r_max
                % Define objective function: Maximize Efficiency
                % Pass 'r' fixed, optimize 'f'
                obj_fun = @(f) -drive_efficiency_circle_helper(r, f, R, T, drive_function);
                
                % --- Loop 2: Continuous optimization for Frequency ---
                % fminbnd efficiently finds the peak frequency for this specific radius
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
        
        % ================== 4. Merge and Save Results ==================
        save_results(results, output_dir, DRIVE_TYPE);
        
    catch ME
        % Error Handling for Command Line
        fprintf('Error occurred:\n%s\n', ME.message);
        fprintf('Stack trace:\n');
        disp(ME.stack);
        exit(1); % Exit with error code
    end
    
    fprintf('Job Complete.\n');
    exit(0); % Exit successfully
end

function save_results(results, output_dir, DRIVE_TYPE)
    output_filename = sprintf('%s/%s_optimal_circle_release.csv', output_dir, DRIVE_TYPE);
    fprintf('\nSaving results to: %s\n', output_filename);
    
    if ~exist(output_dir, 'dir')
       mkdir(output_dir);
    end
    
    % Filter out empty results in case of partial failures
    results = results(~cellfun('isempty', results));
    
    if isempty(results)
        fprintf('No results to save.\n');
        return;
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

function output = drive_efficiency_circle_helper(radius, frequency, R, T, drive_function)
    % FIXED: Use R instead of hardcoded 200 to match the domain size
    intro_vec = zeros(1, R);
    
    % Ensure radius doesn't exceed domain (safety check)
    r_idx = min(round(radius), R);
    
    if r_idx > 0
        intro_vec(1:r_idx) = frequency;
    end
    
    % Run simulation
    output = drive_function(intro_vec, R, T);
end