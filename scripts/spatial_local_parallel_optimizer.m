clear; close all;
%% ================== 1. Configuration ==================
% This script runs MULTIPLE LOCAL optimizations in PARALLEL.
% It reads starting points for different time values from a CSV file.
% --- Define the specific run you want to analyze ---
drive_type = '2-locus_TARE';
% --- Define the TIME points you want to run in parallel ---
T_range = [300 400 500 1000];
% --- Simulation Parameters ---
R = 400; % Total simulation domain radius
r = 200; % Optimization vector dimension (release radius)
%% ================== 2. Setup based on Drive Type ==================
fprintf('Configuring parallel optimizer for drive type: %s\n', drive_type);
switch drive_type
    case 'Wolbachia'
        simulation_func = @Wolbachia_spatial;
        % --- KEY SETTING: Specify the CSV file with the starting vectors ---
        % This file MUST have a 'time' column and 'f1' through 'f200' columns.
        start_points_csv = './start_points/CifAB_start_points.csv';
    case 'TARE'
        simulation_func = @TARE_spatial; 
    case '2-locus_TARE'
        simulation_func = @Two_locus_TARE_spatial;
    case 'TADE_suppression'
        simulation_func = @TADE_suppression_spatial;
    case 'CifAB'
        simulation_func = @CifAB_spatial;     
    otherwise
        error('Drive type "%s" is not recognized.', drive_type);
end
%% ================== 3. Load Starting Points from CSV ==================
if ~exist(start_points_csv, 'file')
    error('Starting points file not found: %s', start_points_csv);
end
fprintf('Loading starting points from: %s\n', start_points_csv);
start_points_table = readtable(start_points_csv);
%% ================== 4. Run Local Optimizations in Parallel ==================
% Check for Parallel Computing Toolbox
if isempty(ver('parallel'))
    error('Parallel Computing Toolbox is not installed.');
end
fprintf('\n--- Starting parallel optimizations for %d time points ---\n', length(T_range));
% Pre-allocate a cell array to store results from workers
results = cell(length(T_range), 1);
parfor i = 1:length(T_range)
    T = T_range(i);
    
    % --- Find and extract the starting vector for the current T ---
    row_idx = find(start_points_table.time == 100, 1);
    if isempty(row_idx)
        fprintf('Warning: No starting point found for T=%d in CSV. Skipping this iteration.\n', T);
        continue; % Skip to the next iteration in the parfor loop
    end
    
    % Convert the table row for f1-f200 to a numeric vector
    start_point_vector = table2array(start_points_table(row_idx, 2:r+1));
    start_point_vector = start_point_vector ./ (1 + start_point_vector);
    % --- Setup and run fmincon on the parallel worker ---
    objective_func = @(x) drive_efficiency(x, R, T, simulation_func);
    lb = zeros(1, r);
    ub = ones(1, r);
    
    % IMPORTANT: Display and Plot functions must be disabled for parfor
    opts = optimoptions(@fmincon, ...
        'Display', 'none', ... % Cannot display iterative output in parfor
        'MaxFunctionEvaluations', 5e5, ...
        'FunctionTolerance', 1e-4, ...
        'StepTolerance', 1e-4);
        
    [x_final, f_final, exitflag] = fmincon(objective_func, start_point_vector, [], [], [], [], lb, ub, [], opts);
    
    % --- Store results in a struct to pass back to the client ---
    % FIX: The struct must be created in a single, atomic assignment
    % to avoid classification errors in a parfor loop.
    results{i} = struct('T', T, 'x_final', x_final, 'f_final', f_final, 'exitflag', exitflag);
    
    % Note: fprintf inside a parfor loop prints to the worker's command window,
    % which may not be visible. It is kept here for potential debugging.
    fprintf('Finished optimization for T = %d on worker.\n', T);
end
fprintf('\n--- All parallel optimizations completed ---\n');
%% ================== 5. Save and Sort All Results ==================
output_filename = sprintf('./results/variable_release/%s_optimal_variable_release.csv', drive_type);
fprintf('\nSaving all results to %s...\n', output_filename);

% --- Ensure the output directory exists ---
output_dir = fileparts(output_filename);
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
   fprintf('Output directory created: %s\n', output_dir);
end

for i = 1:length(results)
    result = results{i};
    
    % Skip if this iteration was skipped (e.g., no start point found)
    if isempty(result)
        continue;
    end
    
    % (The saving logic is the same as your previous script, just looped)
    header = {'time'};
    f_cols = arrayfun(@(k) sprintf('f%d', k), 1:r, 'UniformOutput', false);
    header = [header, f_cols, {'efficiency', 'exitflag'}];
    
    new_data_row = [result.T, result.x_final, -result.f_final, result.exitflag];
    new_data_table = array2table(new_data_row, 'VariableNames', header);
    
    if exist(output_filename, 'file')
        existing_data = readtable(output_filename);
        existing_data(existing_data.time == result.T, :) = [];
        updated_data = [existing_data; new_data_table];
        sorted_data = sortrows(updated_data, 'time');
        writetable(sorted_data, output_filename);
    else
        % If the file doesn't exist, we must handle the first write carefully.
        writetable(new_data_table, output_filename);
    end
end
fprintf('Save complete. CSV file is updated and sorted.\n');
fprintf('======================================================\n');
%% ================== Objective Function Wrapper ==================
function output = drive_efficiency(x, R, T, sim_func)
    % This function will be sent to each parallel worker.
    res = sim_func(x, R, T);
    output = -res;
end

