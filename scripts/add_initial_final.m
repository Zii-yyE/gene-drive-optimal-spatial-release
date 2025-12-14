function add_initial_final(DRIVE_TYPE, RELEASE_TYPE)
%ADD_SIMULATION_STATES Reads an existing optimization results CSV, re-runs
%the simulation for each entry, and appends the 'initial' and 'final'
%state values as two new columns. This script intelligently skips entries
%that already have valid data.
%
%   Usage:
%       add_simulation_states('TARE', 'circle')
%       add_simulation_states('Wolbachia', 'variable')
%
%   Inputs:
%       DRIVE_TYPE (string): The name of the gene drive (e.g., 'TARE').
%       RELEASE_TYPE (string): The type of release pattern ('circle' or 'variable').

%% ================== 1. Configuration ==================
R = 400; % Domain Radius, must match the original simulation
r_dim = 200; % Dimension of the introduction vector

fprintf('--- Starting to add initial/final states for: %s (%s) ---\n', DRIVE_TYPE, RELEASE_TYPE);

%% ================== 2. Setup Paths and Function Handles ==================
% --- Select the appropriate simulation function based on the drive type ---
switch DRIVE_TYPE
    case 'TARE'
        drive_function = @TARE_spatial;
    case '2-locus_TARE'
        drive_function = @Two_locus_TARE_spatial;
    case 'TADE_suppression'
        drive_function = @TADE_suppression_spatial;
    case 'CifAB'
        drive_function = @CifAB_spatial;
    case 'Wolbachia'
        drive_function = @Wolbachia_spatial;
    otherwise
        error("Unknown drive type specified: '%s'", DRIVE_TYPE);
end

% --- Determine input and output file paths based on release type ---
switch RELEASE_TYPE
    case 'circle'
        input_path = sprintf('../results/circle_release/%s_optimal_circle_release.csv', DRIVE_TYPE);
    case 'variable'
        input_path = sprintf('../results/spatial_release/%s_optimal_variable_release.csv', DRIVE_TYPE);
    otherwise
        error("Unknown release type: '%s'. Use 'circle' or 'variable'.", RELEASE_TYPE);
end

% Define a new output path to avoid overwriting original data
output_path = input_path;

%% ================== 3. Load Existing Data & Prepare Columns ==================
if ~exist(input_path, 'file')
    error('Input file not found: %s', input_path);
end
fprintf('Loading data from: %s\n', input_path);
data_table = readtable(input_path);

num_rows = height(data_table);
fprintf('Found %d entries to process.\n', num_rows);

% --- MODIFICATION: Smartly add 'initial' and 'final' columns if they don't exist ---
if ~ismember('initial', data_table.Properties.VariableNames)
    fprintf("Column 'initial' not found. Adding it.\n");
    data_table.initial = nan(num_rows, 1); % Use NaN to mark missing data
end
if ~ismember('final', data_table.Properties.VariableNames)
    fprintf("Column 'final' not found. Adding it.\n");
    data_table.final = nan(num_rows, 1);
end


%% ================== 4. Process Each Entry and Re-run Simulation ==================
for i = 1:num_rows
    % Get the time for the current row
    T = data_table.time(i);
    
    % --- MODIFICATION: Check if data already exists for this row ---
    % A value of 0 might be a valid result, so we also check for NaN (Not-a-Number).
    if ~isnan(data_table.initial(i)) && ~isnan(data_table.final(i)) && data_table.initial(i) ~= 0
        fprintf('Skipping T = %d (%d/%d), data already exists.\n', T, i, num_rows);
        continue; % Skip to the next iteration
    end
    
    fprintf('Processing T = %d (%d/%d)...\n', T, i, num_rows);
    
    % --- Reconstruct the introduction vector based on the release type ---
    intro_vec = zeros(1, r_dim);
    if strcmp(RELEASE_TYPE, 'circle')
        radius = floor(data_table.radius(i));
        frequency = data_table.frequency(i);
        if radius > 0
            intro_vec(1:radius) = frequency;
        end
    else % 'variable'
        % Robustly get f1, f2, ... f200 columns
        var_names = arrayfun(@(k) sprintf('f%d', k), 1:r_dim, 'UniformOutput', false);
        intro_vec = data_table{i, var_names};
    end
    
    % --- Re-run the simulation ---
    % FIX: You must capture all three outputs from the function into separate
    % variables. The tilde (~) is used to ignore the first output (res),
    % as this script does not need it.
    [~, initial_val, final_val] = drive_function(intro_vec, R, T);
    
    % --- Store results directly into the table ---
    data_table.initial(i) = initial_val;
    data_table.final(i) = final_val;
end

%% ================== 5. Save Updated Table ==================
fprintf('\nAll simulations complete. Saving updated table...\n');

% Save the updated table to the new output file
writetable(data_table, output_path);

fprintf('Successfully saved enriched data to: %s\n', output_path);
fprintf('======================================================\n');

end

