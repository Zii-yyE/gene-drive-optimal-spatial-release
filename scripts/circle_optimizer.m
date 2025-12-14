function circle_optimizer(DRIVE_TYPE)
% This script performs a grid search to find the optimal circle release
% parameters (radius and frequency) for a given DRIVE_TYPE.
% It runs the grid search for multiple simulation times in PARALLEL and
% intelligently merges the results into a single, sorted CSV file.

%% ================== 1. Main Configuration ==================
R = 400; % Domain Radius
r_max = 200; % Max Introduction Radius
num_radius_steps = 101; % Number of radius points to test
num_frequency_steps = 100; % Number of frequency points to test

output_dir = './raw_results/circle_release'; % Define the output directory

fprintf('--- Starting Parallel Grid Search for: %s ---\n', DRIVE_TYPE);

%% ================== 2. Select Drive Function & Time Range ==================
% Get the correct function handle and time range based on the drive type.
switch DRIVE_TYPE
    case 'TARE'
        drive_function = @TARE_spatial;
        T_range = 10:10:100;
    case 'Two_Locus_TARE'
        drive_function = @Two_locus_TARE_spatial;
        T_range = 10:10:100;
    case 'TADE_suppression'
        drive_function = @TADE_suppression_spatial;
        T_range = [10 50 75 100 125 150 175 200];
    case 'CifAB'
        drive_function = @CifAB_spatial;
        T_range = [10 20 30 40 50];
    case 'Wolbachia'
        drive_function = @Wolbachia_spatial;
        T_range = 10:10:100;
    otherwise
        error("Unknown drive type specified: '%s'", DRIVE_TYPE);
end

%% ================== 3. Run Grid Searches in Parallel for Each Time ==================
fprintf('\n--- Processing %d time points in parallel ---\n', length(T_range));

% Pre-allocate a cell array to store results from workers
results = cell(length(T_range), 1);
tic; % Start a timer for the entire parallel operation

% --- MODIFICATION: Changed to a parfor loop for parallel execution ---
parfor t_idx = 1:length(T_range)
    T = T_range(t_idx);
    
    % Define the grid for this specific worker
    radius_grid = linspace(0, r_max, num_radius_steps);
    frequency_grid = linspace(0, 0.99, num_frequency_steps);
    
    % Pre-allocate the landscape matrix for this worker
    landscape = zeros(num_radius_steps, num_frequency_steps);
    
    % This inner loop runs on the parallel worker
    for i = 1:num_radius_steps
        temp_slice = zeros(1, num_frequency_steps);
        for j = 1:num_frequency_steps
            params = [radius_grid(i), frequency_grid(j)];
            temp_slice(j) = drive_efficiency_circle(params, R, T, drive_function);         
        end
        landscape(i, :) = temp_slice;
    end
    
    % --- Find the best result for this T ---
    [max_efficiency, linear_idx] = max(landscape(:));
    [i_max, j_max] = ind2sub(size(landscape), linear_idx);
    
    optimal_radius = radius_grid(i_max);
    optimal_frequency = frequency_grid(j_max);
    
    % --- Store results in a struct ---
    results{t_idx} = struct('time', T, 'radius', optimal_radius, 'frequency', optimal_frequency, 'efficiency', max_efficiency);
    fprintf('Finished grid search for T = %d on worker.\n', T);
end
total_elapsed_time = toc;
fprintf('\n--- All parallel grid searches completed in %.2f seconds ---\n', total_elapsed_time);

%% ================== 4. Merge and Save All Results to a Single CSV ==================
output_filename = sprintf('%s/%s_optimal_circle_release.csv', output_dir, DRIVE_TYPE);
fprintf('\nSaving/updating results in: %s\n', output_filename);

% --- Ensure the output directory exists ---
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
   fprintf('Output directory created: %s\n', output_dir);
end

% --- Consolidate new results into a table ---
new_results_table = struct2table(cat(1, results{:}));

% --- MODIFICATION: Robust saving logic ---
if exist(output_filename, 'file')
    % If file exists, read it, merge new data, and re-sort
    fprintf('File exists. Reading, merging, and re-sorting...\n');
    existing_data = readtable(output_filename);
    
    % Get the times from the new results
    new_times = new_results_table.time;
    
    % Remove any old rows that have the same times as the new results
    existing_data(ismember(existing_data.time, new_times), :) = [];
    
    % Combine the old (filtered) data with the new data
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

function output = drive_efficiency_circle(params, R, T, drive_function)
% Helper function to prepare the input vector and call the simulation.
radius = floor(params(1)); % Ensure radius is an integer.
frequency = params(2);

% Assuming a fixed optimization vector size of 200
intro_vec = zeros(1, 200); 
if radius > 0
    intro_vec(1:radius) = frequency;
end

% Call the appropriate drive function passed as a handle.
output = drive_function(intro_vec, R, T);
end
