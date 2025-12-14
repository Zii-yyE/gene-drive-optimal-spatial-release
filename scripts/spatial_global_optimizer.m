clear; close all;
%% ================== 1. Configuration ==================
drive_type = 'CifAB';
T_range = [100];
R = 400; % Total simulation domain radius
r = 200; % Optimization vector dimension (release radius)

%% ================== 2. Setup based on Drive Type ==================
fprintf('Configuring optimizer for drive type: %s\n', drive_type);

custom_start_point = []; % This will hold an optional, manually defined start point

switch drive_type
    case 'Wolbachia'
        simulation_func = @Wolbachia_spatial;
        optimal_circle_csv = '../results/circle_release/Wolbachia_optimal_circle_release.csv';

        custom_start_point = zeros(1, r);
        custom_start_point(1:11) = [0.27712081503271 , 0.274911463699575, 0.270707808309674, 0.263779919826942, 0.25370402694546 , 0.240128345518471, 0.222571399186681, 0.200581917034892 , 0.173554624818741  , 0.130004879723616   , 0.0266426595140543];

    case 'TARE'
        simulation_func = @TARE_spatial; 
        optimal_circle_csv = '../results/circle_release/TARE_optimal_circle_release.csv';

    % case '2-locus_TARE'
    %     simulation_func = @2-locus_TARE_spatial;
    %     optimal_circle_csv = '../results/circle_release/2-locus_TARE_optimal_circle_release.csv';

    case 'TADE_suppression'
        simulation_func = @TADE_suppression_spatial;
        optimal_circle_csv = '../results/circle_release/TADE_suppression_optimal_circle_release.csv';
        
    otherwise
        error('Drive type "%s" is not recognized. Please add a new case to the switch statement in section 2.', drive_type);
end

%% ================== 3. Main Optimization Loop ==================
% Loop over each simulation time T to find the optimal release pattern.
for i = 1:length(T_range)
    T = T_range(i);
    fprintf('\n--- Starting optimization for T = %d ---\n', T);
    %% ================== 4. Start Point Generation ==================
    % --- 4.1. Generate random starting points ---
    n_random_pts = 1;
    random_pts_matrix = zeros(n_random_pts, r);
    for j = 1:n_random_pts
        random_values = rand(1, 10);
        random_pts_matrix(j, :) = repelem(random_values, 20);
    end

    all_points = random_pts_matrix; % Start with random points

    % --- 4.2. Add the optimal circle release as a smart starting point ---
    if exist(optimal_circle_csv, 'file')
        data = readtable(optimal_circle_csv);
        optimal_row = data(data.time == T, :);
        if ~isempty(optimal_row)
            radius = optimal_row.radius(1);
            frequency = optimal_row.frequency(1);
            optimal_circle_point = zeros(1, r);
            optimal_circle_point(1:radius) = frequency;
            all_points = [optimal_circle_point; all_points];
            fprintf('Successfully added optimal circle release as a starting point.\n');
        else
            fprintf('Warning: T=%d not found in %s. Using random points only.\n', T, optimal_circle_csv);
        end
    else
        fprintf('Warning: %s not found. Using random points only.\n', optimal_circle_csv);
    end
    
    % ----- NEW FEATURE: Add the custom start point if it was defined -----
    if ~isempty(custom_start_point)
        all_points = [all_points; custom_start_point]; % Append the custom point
        fprintf('Successfully added manually defined custom start point.\n');
    end
    
    tpoints = CustomStartPointSet(all_points);
    %% ================== 5. Optimization Problem Setup ==================
    % --- 5.1. MultiStart Global Solver Settings ---
    ms = MultiStart;
    ms.UseParallel = true;
    ms.Display = 'iter';
    ms.FunctionTolerance = 1e-3;
    ms.XTolerance = 1e-3;
    ms.MaxTime = 500000;
    ms.StartPointsToRun = 'bounds';
    % --- 5.2. fmincon Local Solver Settings ---
    lb = zeros(1, r);
    ub = ones(1, r);
    opts = optimoptions(@fmincon, "Display", "none");
    opts.MaxFunctionEvaluations = 5e5;
    % --- 5.3. Problem Definition ---
    objective_func = @(x) drive_efficiency(x, R, T, simulation_func);
    x0 = all_points(1,:);
    problem = createOptimProblem('fmincon', 'x0', x0, ...
        'objective', objective_func, 'lb', lb, 'ub', ub, 'options', opts);
    %% ================== 6. Run Optimization and Save Results ==================
    rng(1); % for reproducibility
    [xming, fming, flagg] = run(ms, problem, tpoints);
    % --- Save the best result ---
    file_path = sprintf('../results/spatial_release/%s_R%d_r%d_T%d.csv', drive_type, R, r, T);
    
    output_dir = fileparts(file_path);
    if ~exist(output_dir, 'dir')
       mkdir(output_dir);
       fprintf('Output directory created: %s\n', output_dir);
    end
    
    if ~exist(file_path, 'file')
        header = {'time'};
        f_cols = arrayfun(@(i) sprintf('f%d', i), 1:r, 'UniformOutput', false);
        header = [header, f_cols, {'efficiency', 'exitflag'}];
        writecell(header, file_path);
    end
    data_row = [T, xming, -fming, flagg];
    writematrix(data_row, file_path, "WriteMode", "append");
    fprintf('Result for T=%d saved to %s\n', T, file_path);
    
end
fprintf('\nAll optimizations completed.\n');
%% ================== Objective Function Wrapper ==================
function output = drive_efficiency(x, R, T, sim_func)
    res = sim_func(x, R, T);
    output = -res;
end

