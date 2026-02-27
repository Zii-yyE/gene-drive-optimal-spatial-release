function run_TARE_optimization()
    try
        % ================== 1. Configuration ==================
        drive_type = 'TARE';
        T_range = [10000];
        R = 10000; 
        r = 200; 
        start_points_csv = './start_points/TARE_start_points.csv';
        
        % ================== 2. Setup ==================
        fprintf('Configuring optimizer for: %s\n', drive_type);
        simulation_func = @TARE_spatial;
        
        if ~exist(start_points_csv, 'file')
            error('Starting points file not found: %s', start_points_csv);
        end
        start_points_table = readtable(start_points_csv);
        
        % ================== 3. Parallel Pool Setup ==================
        delete(gcp('nocreate')); 
        parpool(); 
        
        fprintf('\n--- Starting parallel optimizations ---\n');
        results = cell(length(T_range), 1);
        
        % ================== 4. Run Optimization ==================
        for i = 1:length(T_range)
            maxNumCompThreads(1); 
            
            T = T_range(i);
            
            row_idx = find(start_points_table.time == 10000, 1);
            
            if isempty(row_idx)
                fprintf('Skipping T=%d (No start point)\n', T);
                continue;
            end
            
            start_vec = table2array(start_points_table(row_idx, 2:r+1));
            start_vec = start_vec ./ (1 + start_vec);
            
            objective_func = @(x) drive_efficiency_wrapper(x, R, T, simulation_func);
            lb = zeros(1, r);
            ub = ones(1, r);
            
            opts = optimoptions(@fmincon, ...
                'Display', 'none', ...
                'MaxFunctionEvaluations', 5e5, ...
                'FunctionTolerance', 1e-4, ...
                'StepTolerance', 1e-4, ...
                'Algorithm', 'sqp');
                
            fprintf('Worker started optimizing T=%d...\n', T);
            [x_final, f_final, exitflag] = fmincon(objective_func, start_vec, [], [], [], [], lb, ub, [], opts);
            
            results{i} = struct('T', T, 'x_final', x_final, 'f_final', f_final, 'exitflag', exitflag);
            fprintf('Worker finished T=%d.\n', T);
        end
        
        % ================== 5. Save Results ==================
        save_results(results, drive_type, r);
        
    catch ME
        fprintf('Error occurred:\n%s\n', ME.message);
        fprintf('Stack trace:\n');
        disp(ME.stack);
        exit(1); 
    end
    
    fprintf('Job Complete.\n');
end

function output = drive_efficiency_wrapper(x, R, T, sim_func)
    output = -sim_func(x, R, T);
end

function save_results(results, drive_type, r)
    output_filename = sprintf('./results/variable_release/%s_optimal_variable_release.csv', drive_type);
    output_dir = fileparts(output_filename);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    for i = 1:length(results)
        result = results{i};
        if isempty(result), continue; end
        
        header = {'time'};
        f_cols = arrayfun(@(k) sprintf('f%d', k), 1:r, 'UniformOutput', false);
        header = [header, f_cols, {'efficiency', 'exitflag'}];
        
        new_row = [result.T, result.x_final, -result.f_final, result.exitflag];
        new_table = array2table(new_row, 'VariableNames', header);
        
        if exist(output_filename, 'file')
            existing = readtable(output_filename);
            existing(existing.time == result.T, :) = [];
            final = sortrows([existing; new_table], 'time');
            writetable(final, output_filename);
        else
            writetable(new_table, output_filename);
        end
    end
    fprintf('Results saved to %s\n', output_filename);
end
