function add_initial_final_releasesize(DRIVE_TYPE, RELEASE_TYPE)
%Reads an existing optimization results CSV, re-runs the simulation 
% to append 'initial', 'final', and 'release_size' columns.
%
%   Usage:
%       add_initial_final_releasesize('TARE', 'circle')
%       add_initial_final_releasesize('Wolbachia', 'variable')

%% ================== 1. Configuration ==================
R = 400; % Domain Radius
r_dim = 200; % Dimension of intro vector
fprintf('--- Updating data for: %s (%s) ---\n', DRIVE_TYPE, RELEASE_TYPE);

%% ================== 2. Setup Paths ==================
switch DRIVE_TYPE
    case 'TARE', drive_function = @TARE_spatial;
    case '2-locus_TARE', drive_function = @Two_locus_TARE_spatial;
    case 'TADE_suppression', drive_function = @TADE_suppression_spatial;
    case 'CifAB', drive_function = @CifAB_spatial;
    case 'Wolbachia', drive_function = @Wolbachia_spatial;
    otherwise, error("Unknown drive type: '%s'", DRIVE_TYPE);
end

switch RELEASE_TYPE
    case 'circle'
        input_path = sprintf('../results/circle_release/%s_optimal_circle_release.csv', DRIVE_TYPE);
    case 'variable'
        input_path = sprintf('../results/variable_release/%s_optimal_variable_release.csv', DRIVE_TYPE);
    otherwise, error("Unknown release type: '%s'", RELEASE_TYPE);
end

output_path = input_path; % Overwrite mode

%% ================== 3. Load & Initialize Columns ==================
if ~exist(input_path, 'file'), error('File not found: %s', input_path); end

data_table = readtable(input_path);
num_rows = height(data_table);
fprintf('Loaded %d entries.\n', num_rows);

% --- Initialize columns if they don't exist (Fix for Issue 3) ---
if ~ismember('initial', data_table.Properties.VariableNames)
    data_table.initial = nan(num_rows, 1);
end
if ~ismember('final', data_table.Properties.VariableNames)
    data_table.final = nan(num_rows, 1);
end
% Use 'release_size' instead of 'size' to avoid conflict with built-in function
if ~ismember('release_size', data_table.Properties.VariableNames)
    fprintf("Column 'release_size' added.\n");
    data_table.release_size = nan(num_rows, 1);
end

%% ================== 4. Process Loop ==================
for i = 1:num_rows
    T = data_table.time(i);
    
    % --- Fix for Issue 1: Check if release_size is missing ---
    % Only skip if ALL data is present. If size is NaN, we must run.
    has_initial = ~isnan(data_table.initial(i)) && data_table.initial(i) ~= 0;
    has_final   = ~isnan(data_table.final(i));
    has_size    = ~isnan(data_table.release_size(i));
    
    if has_initial && has_final && has_size
        % verify: print less frequently to reduce clutter
        % fprintf('Skipping T = %d (Data complete).\n', T);
        continue; 
    end
    
    fprintf('Processing T = %d (%d/%d)...\n', T, i, num_rows);
    
    % --- Reconstruct Intro Vector ---
    intro_vec = zeros(1, r_dim);
    if strcmp(RELEASE_TYPE, 'circle')
        radius = floor(data_table.radius(i));
        frequency = data_table.frequency(i);
        if radius > 0
            % FIX: Ensure radius doesn't exceed dimension
            r_idx = min(radius, r_dim);
            intro_vec(1:r_idx) = frequency;
        end
    else % variable
        var_names = arrayfun(@(k) sprintf('f%d', k), 1:r_dim, 'UniformOutput', false);
        intro_vec = data_table{i, var_names};
    end
    
    % --- Re-run Simulation ---
    % Correction for Issue 2: Use 4 output arguments as originally intended.
    % Assuming your spatial functions return: [res, initial_val, final_val, release_size]
    [~, initial_val, final_val, release_size] = drive_function(intro_vec, R, T);
    
    % --- Store Results ---
    data_table.initial(i) = initial_val;
    data_table.final(i)   = final_val;
    data_table.release_size(i) = release_size; 
end

%% ================== 5. Save ==================
writetable(data_table, output_path);
fprintf('Success. Updated file: %s\n', output_path);
end