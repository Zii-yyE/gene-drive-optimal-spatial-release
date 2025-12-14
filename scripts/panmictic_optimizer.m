function panmictic_optimizer(DRIVE_TYPE)
% This script finds the optimal release frequency for a panmictic model,
% which is equivalent to a spatial model with a fixed, global release.
% It runs the search for multiple simulation times in PARALLEL and
% intelligently merges the results into a single, sorted CSV file.

%% ================== 1. Main Configuration ==================
R = 400; % Domain Radius
% MODIFICATION: Removed r_max. For a panmictic model, the release radius
% is always equal to the full domain radius, R.
num_frequency_steps = 200; % Number of frequency points to test

output_dir = '../results/panmictic_release'; % Define the output directory

fprintf('--- Starting Panmictic Optimization for: %s ---\n', DRIVE_TYPE);

%% ================== 2. Select Drive Function & Time Range ==================
% Get the correct function handle and time range based on the drive type.
switch DRIVE_TYPE
    case 'TARE'
        drive_function = @TARE_spatial;
        T_range = 10:10:100;
    case '2-locus_TARE'
        drive_function = @Two_locus_TARE_spatial;
        T_range = 10:10:100;
    case 'TADE_suppression'
        drive_function = @TADE_suppression_spatial;
        T_range = [10 50 75 100 125 150 175 200];
    case 'CifAB'
        drive_function = @CifAB_spatial;
        T_range = [1000];
    case 'Wolbachia'
        drive_function = @Wolbachia_spatial;
        T_range = 10:10:100;
    otherwise
        error("Unknown drive type specified: '%s'", DRIVE_TYPE);
end

%% ================== 3. Run Searches in Parallel for Each Time ==================
fprintf('\n--- Processing %d time points in parallel ---\n', length(T_range));

% Pre-allocate a cell array to store results from workers
results = cell(length(T_range), 1);
tic; % Start a timer for the entire parallel operation

parfor t_idx = 1:length(T_range)
    T = T_range(t_idx);
    
    % Define the 1D grid for this specific worker
    frequency_grid = linspace(0, 0.99, num_frequency_steps);
    
    % Pre-allocate an efficiency vector for this worker
    efficiency_results = zeros(1, num_frequency_steps);
    
    % This inner loop runs on the parallel worker
    for j = 1:num_frequency_steps
        frequency = frequency_grid(j);
        % MODIFICATION: Removed r_max from the function call.
        efficiency_results(j) = drive_efficiency_panmictic(frequency, R, T, drive_function);
    end
    
    % --- Find the best result for this T ---
    [max_efficiency, j_max] = max(efficiency_results);
    optimal_frequency = frequency_grid(j_max);
    
    % --- Store results in a struct ---
    results{t_idx} = struct('time', T, 'frequency', optimal_frequency, 'efficiency', max_efficiency);
    fprintf('Finished panmictic search for T = %d on worker.\n', T);
end
total_elapsed_time = toc;
fprintf('\n--- All parallel searches completed in %.2f seconds ---\n', total_elapsed_time);

%% ================== 4. Merge and Save All Results to a Single CSV ==================
output_filename = sprintf('%s/%s_optimal_panmictic_release.csv', output_dir, DRIVE_TYPE);
fprintf('\nSaving/updating results in: %s\n', output_filename);

% --- Ensure the output directory exists ---
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
   fprintf('Output directory created: %s\n', output_dir);
end

% --- Consolidate new results into a table ---
new_results_table = struct2table(cat(1, results{:}));

% --- Robust saving logic ---
if exist(output_filename, 'file')
    % If file exists, read it, merge new data, and re-sort
    fprintf('File exists. Reading, merging, and re-sorting...\n');
    existing_data = readtable(output_filename);
    
    new_times = new_results_table.time;
    existing_data(ismember(existing_data.time, new_times), :) = [];
    final_data_to_save = [existing_data; new_results_table];
else
    % If file does not exist, the final data is just the new results
    fprintf('File does not exist. Creating new file...\n');
    final_data_to_save = new_results_table;
end

% --- Sort the final, combined table by time and write to file ---
sorted_table = sortrows(final_data_to_save, 'time');
writetable(sorted_table, output_filename);

fprintf('Save complete. CSV is updated and sorted.\n');
fprintf('======================================================\n');
end

function output = drive_efficiency_panmictic(frequency, R, T, drive_function)
% Helper function for the panmictic model simulation.
% MODIFICATION: The radius is now fixed to the full domain size R.
% The intro_vec size is now determined by R, not a hardcoded value.
intro_vec = zeros(1, R);
intro_vec(1:R) = frequency;

% Call the appropriate drive function passed as a handle.
output = drive_function(intro_vec, R, T);
end

