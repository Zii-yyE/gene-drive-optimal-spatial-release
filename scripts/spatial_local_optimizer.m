clear; close all;
%% ================== 1. Configuration ==================
% This script runs a SINGLE LOCAL optimization (fmincon) from a
% manually defined starting point and provides detailed output.

% --- Define the specific run you want to analyze ---
drive_type = 'Wolbachia';
T = 20; % A single time value, not a range

% --- KEY SETTING: Define your 200-element starting vector here ---
% This is where you provide your custom initial release pattern.
% It must be a 1-by-200 row vector.
% Example 1: A "center" release
user_defined_start_point = 1e-8*ones(1, 200);
user_defined_start_point(1:200) = [0.598770374022481,0.484664988742014,0.454249659820879,0.473206232483993,0.49345091015734,0.493591873101838,0.494584591308686,0.500059825451692,0.503936171368765,0.512945628440681,0.517160817035796,0.519944776801473,0.526017607258081,0.531439625515603,0.536708375839986,0.545603125078256,0.554227703452948,0.56521090268916,0.577788710462723,0.593166404015314,0.612459076435917,0.632520845962793,0.657404574729883,0.685767088626093,0.719216654084475,0.749523115334616,0.785669802048172,0.824059205066584,0.862807967921815,0.900617315008194,0.93777249128937,0.975810361300103,1.00196169363739,1.02794717010933,1.03948248269369,1.03286641601526,1.00243409751905,0.932052429122129,0.781627304146769,0.425505369305683,0.00267745924679446,3.8711541563599e-4,1.7444098807196e-4,1.0758679028777e-4,7.78693697922459e-5,6.20397830585614e-5,5.26622482410776e-5,4.67152774401493e-5,4.278524246933e-5,4.01413475934815e-5,3.83915156776471e-5,3.73196466995233e-5,3.68278208001786e-5,3.69067079406416e-5,3.76334927461499e-5,3.91905057388966e-5,4.19282910618966e-5,4.65105793199227e-5,5.42686205306139e-5,6.81529160496825e-5,9.59818197413596e-5,1.6514954535425e-4,4.4357483773176e-4,0.0849053374751395,0.555489824516137,0.767314697645265,0.857570019668455,0.900782842156152,0.914688042638734,0.910393327275628,0.892658686650791,0.871600814292153,0.844317644454471,0.82047183071447,0.798465467761268,0.782970434927712,0.775124954398461,0.774287951718631,0.781822986420327,0.797644262140001,0.82052198423333,0.848049768265849,0.879243520475697,0.912759899649938,0.943789907306771,0.969241464976355,0.984824909917325,0.986059533890499,0.962543142359813,0.904125472773708,0.775904090997989,0.46359261896878,0.00199423532573948,1.9261320750438e-4,8.38270772830764e-5,5.12429115695123e-5,3.70423628612258e-5,2.9571633953458e-5,2.51843649413839e-5,2.24226880752854e-5,2.06059186320961e-5,1.9381582341559e-5,1.85522034745864e-5,1.80059361702972e-5,1.76788963601485e-5,1.7542541173349e-5,1.75933992206714e-5,1.78552516885765e-5,1.83824033940453e-5,1.92731058912461e-5,2.07021415124744e-5,2.29838403978234e-5,2.67279930510203e-5,3.32580631980854e-5,4.59487238514983e-5,7.59654308478509e-5,1.8352567484519e-4,0.00671297255318484,0.503035817739599,0.758625707886021,0.864719883955022,0.915372843004334,0.933450275064736,0.931116524325151,0.915177831476594,0.891789144188403,0.864461387423745,0.837482016772803,0.812841464721162,0.793733355831844,0.78105534579867,0.776343343974411,0.780391083046585,0.792028446743778,0.811713795853763,0.837247064037081,0.866813202371717,0.897837383650969,0.928758220678453,0.955100629591718,0.972706110905869,0.976104743661043,0.957223416166995,0.903895424334102,0.794071539289066,0.530692388176736,0.00927973029730025,1.4403041885648e-4,5.81189871296231e-5,3.45539654109278e-5,2.46405938945324e-5,1.95278599691766e-5,1.65645557446587e-5,1.4718569989385e-5,1.35153617828077e-5,1.27111225852189e-5,1.21706076231975e-5,1.18161779480136e-5,1.16028911425158e-5,1.15083431326055e-5,1.15268469941601e-5,1.16677056224213e-5,1.19582652679951e-5,1.2449197872123e-5,1.32284488238779e-5,1.44511686056937e-5,1.64047258928553e-5,1.9671289561452e-5,2.5597375415927e-5,3.79654827665541e-5,7.14268700012343e-5,2.4472143633229e-4,0.221380141893841,0.645081482027434,0.822230188125206,0.898296339982147,0.932359410022497,0.939554949489343,0.93015134033998,0.909221310605832,0.882869476046514,0.854223891636587,0.826776122320635,0.803575304428646,0.786370868816542,0.775927039892999,0.773438205035254,0.780361984234859,0.795498956272692,0.817614779873778,0.845426653201938,0.876721154062993,0.90930815904618,0.940086479523551,0.965691909744389,0.980432937834024,0.977938911673628,0.94443179129806,0.902924713707742, 1.027840141];
user_defined_start_point = user_defined_start_point ./ (1+user_defined_start_point);
% --- Simulation Parameters (should match your other script) ---
R = 400; % Total simulation domain radius
r = 200; % Optimization vector dimension (release radius)


%% ================== 2. Setup based on Drive Type ==================
% Select the appropriate simulation function based on the drive_type.
fprintf('Configuring local optimizer for drive type: %s\n', drive_type);
switch drive_type
    case 'Wolbachia'
        simulation_func = @Wolbachia_spatial;
    case 'TARE'
        simulation_func = @TARE_spatial;
    case '2-locus_TARE'
        simulation_func = @Two-locus_TARE_spatial;
    case 'TADE_suppression'
        simulation_func = @TADE_suppression_spatial;
    case 'CifAB'
        simulation_func = @CifAB_spatial;
    otherwise
        error('Drive type "%s" is not recognized.', drive_type);
end


%% ================== 3. Optimization Problem Setup ==================
% The objective function is a wrapper that calls the correct simulation.
objective_func = @(x) drive_efficiency(x, R, T, simulation_func);

% Define bounds for the release pattern (frequency from 0 to 1).
lb = zeros(1, r);
ub = ones(1, r);

% --- KEY FEATURE: Configure fmincon for DETAILED iterative output ---
% 'Display', 'iter-detailed': Shows detailed progress for each iteration.
% 'PlotFcn', ... : Creates real-time plots showing the optimization's convergence.
opts = optimoptions(@fmincon, ...
    'Display', 'iter-detailed', ...
    'PlotFcn', {@optimplotfval, @optimplotstepsize, @optimplotfirstorderopt}, ...
    'MaxFunctionEvaluations', 5e5, ...
    'FunctionTolerance', 1e-4, ...
    'StepTolerance', 1e-4);


%% ================== 3.5 Evaluate Initial Starting Point ==================
% ----- NEW FEATURE: Calculate and display the efficiency of the starting point BEFORE optimization -----
fprintf('\n------------------------------------------------------\n');
fprintf('       EVALUATING INITIAL STARTING POINT\n');
fprintf('------------------------------------------------------\n');
initial_fval = objective_func(user_defined_start_point);
% ----- MODIFICATION: Format output to 6 decimal places -----
fprintf('Initial Efficiency of Starting Point: %.6f\n', -initial_fval);
fprintf('------------------------------------------------------\n');


%% ================== 4. Run Local Optimization ==================
fprintf('\n--- Starting local optimization for T = %d ---\n', T);
fprintf('Using the manually defined starting point.\n\n');

% Call fmincon directly to start the optimization.
[x_final, f_final, exitflag, output] = fmincon(objective_func, user_defined_start_point, [], [], [], [], lb, ub, [], opts);


%% ================== 5. Display Final Results Summary ==================
fprintf('\n\n======================================================\n');
fprintf('           LOCAL OPTIMIZATION SUMMARY\n');
fprintf('======================================================\n');

fprintf('Drive Type:                  %s\n', drive_type);
fprintf('Simulation Time (T):         %d\n', T);
fprintf('------------------------------------------------------\n');

% Display the detailed output structure from fmincon
disp(output);

fprintf('Exit Flag:                   %d\n', exitflag);
% Provide an interpretation of the exit flag
if exitflag > 0
    fprintf('   (Meaning: Optimization converged successfully.)\n');
elseif exitflag == 0
    fprintf('   (Meaning: Stopped by iteration or function evaluation limit.)\n');
else
    fprintf('   (Meaning: Optimization did not converge.)\n');
end

% ----- MODIFICATION: Format output to 6 decimal places -----
fprintf('Initial Efficiency:          %.6f\n', -initial_fval);
fprintf('Final Optimal Efficiency:      %.6f\n', -f_final);
fprintf('------------------------------------------------------\n');
fprintf('The initial vector (user_defined_start_point) and the final\n');
fprintf('optimized vector (x_final) are available in the workspace for plotting.\n');
fprintf('======================================================\n');


%% ================== 6. Save and Sort Results to CSV ==================
% ----- NEW FEATURE: Save the result to a CSV, creating or updating it as needed -----
output_filename = sprintf('%s_optimal_variable_release.csv', drive_type);
fprintf('\nSaving result to %s...\n', output_filename);

% Create the header for the CSV file
header = {'time'};
f_cols = arrayfun(@(i) sprintf('f%d', i), 1:r, 'UniformOutput', false);
header = [header, f_cols, {'efficiency', 'exitflag'}];

% Create a table for the new result row
new_data_row = [T, x_final, -f_final, exitflag];
new_data_table = array2table(new_data_row, 'VariableNames', header);

if exist(output_filename, 'file')
    % If file exists, read it, update it, sort it, and write it back
    fprintf('File exists. Reading, updating, and re-sorting...\n');
    existing_data = readtable(output_filename);
    
    % Remove any old row that has the same time T
    existing_data(existing_data.time == T, :) = [];
    
    % Append the new data
    updated_data = [existing_data; new_data_table];
    
    % Sort the combined data by the 'time' column
    sorted_data = sortrows(updated_data, 'time');
    
    % Write the final sorted table to the file
    writetable(sorted_data, output_filename);
else
    % If file does not exist, just write the new data (with header)
    fprintf('File does not exist. Creating new file...\n');
    writetable(new_data_table, output_filename);
end

fprintf('Save complete. CSV file is updated and sorted.\n');
fprintf('======================================================\n');


%% ================== Objective Function Wrapper ==================
function output = drive_efficiency(x, R, T, sim_func)
    % This generic wrapper calls the specified simulation function handle.
    res = sim_func(x, R, T);
    
    % Negate the result because fmincon performs minimization, 
    % but we want to maximize the efficiency.
    output = -res;
end

