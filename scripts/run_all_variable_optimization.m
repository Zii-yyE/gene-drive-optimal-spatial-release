function run_all_variable_optimization()
    % All-in-One Optimization Script for 5 Gene Drives (Variable Release)
    % Runs multiple drive types in PARALLEL for T=10000.
    % MODIFICATIONS:
    % - Optimization dimension (r) reduced to 50 (Truncates start points).
    % - fmincon tolerance reduced to 1e-3 for speed.
    
    try
        % ================== 1. 任务配置表 (自定义区域) ==================
        % 格式: {DriveType, SeedTime, TargetTime, StartPointsFile}
        % SeedTime: 从 CSV 中找哪一行作为起点 (比如用 T=100 的结果跑 T=10000)
        jobs_config = {
            'TARE',             10000,   10000, './start_points/TARE_start_points.csv';
            '2-locus_TARE',     1000,   10000, './start_points/2-locus_TARE_start_points.csv';
            'TADE_suppression', 10000,   10000, './start_points/TADE_suppression_start_points.csv';
            'CifAB',            10000,   10000, './start_points/CifAB_start_points.csv';
            'Wolbachia',        10000,   10000, './start_points/Wolbachia_start_points.csv';
        };
        
        % 全局参数
        R = 10000;  % Domain Radius (防止波跑到边界)
        r = 50;     % [MODIFIED] Optimization vector size reduced to 50
        
        fprintf('--- Starting All-in-One Variable Optimization for %d Drives (T=%d, r=%d) ---\n', ...
                size(jobs_config, 1), 10000, r);

        % ================== 2. 并行池设置 ==================
        delete(gcp('nocreate')); 
        num_workers = size(jobs_config, 1);
        parpool(num_workers);
        
        fprintf('\n--- Parallel pool started. Launching workers... ---\n');
        
        % ================== 3. 并行执行所有任务 ==================
        parfor i = 1:size(jobs_config, 1)
            % [关键] 强制单线程，防止计算资源争抢
            maxNumCompThreads(1);
            
            % --- 解析当前任务参数 ---
            drive_type = jobs_config{i, 1};
            seed_T     = jobs_config{i, 2};
            target_T   = jobs_config{i, 3};
            csv_path   = jobs_config{i, 4};
            
            fprintf('[Worker %d] Preparing %s (Target T=%d)...\n', i, drive_type, target_T);
            
            % --- 1. 选择仿真函数 ---
            sim_func = get_drive_function(drive_type);
            
            % --- 2. 读取起始点 ---
            if ~exist(csv_path, 'file')
                fprintf('[Error] %s: Start CSV not found: %s\n', drive_type, csv_path);
                continue; 
            end
            
            start_table = readtable(csv_path);
            
            % 寻找种子时间点
            row_idx = find(start_table.time == seed_T, 1);
            
            if isempty(row_idx)
                fprintf('[Error] %s: Seed T=%d not found in CSV.\n', drive_type, seed_T);
                continue;
            end
            
            start_vec = table2array(start_table(row_idx, 2:r+1))+0.01;
                        
            % --- 3. 配置 fmincon ---
            obj_fun = @(x) -sim_func(x, R, target_T); 
            lb = zeros(1, r);
            ub = ones(1, r);
            
            % [MODIFIED] 精度调整为 1e-3
            opts = optimoptions(@fmincon, ...
                'Display', 'none', ... 
                'MaxFunctionEvaluations', 5e5, ...
                'FunctionTolerance', 1e-3, ... % Reduced tolerance
                'StepTolerance', 1e-3, ...     % Reduced tolerance
                'Algorithm', 'sqp');
            
            % --- 4. 开始优化 ---
            fprintf('[Worker %d] Running variable optimization for %s...\n', i, drive_type);
            
            t_start = tic;
            [x_final, f_final, exitflag] = fmincon(obj_fun, start_vec, [], [], [], [], lb, ub, [], opts);
            duration = toc(t_start);
            
            % --- 5. 结果保存 ---
            res_struct = struct('T', target_T, 'x_final', x_final, 'f_final', f_final, 'exitflag', exitflag);
            
            save_single_result(res_struct, drive_type, r);
            
            fprintf('[Worker %d] Finished %s in %.1fs. Efficiency: %.4f\n', ...
                    i, drive_type, duration, -f_final);
        end
        
    catch ME
        fprintf('Global Error occurred:\n%s\n', ME.message);
        disp(ME.stack);
        exit(1);
    end
    
    fprintf('All jobs complete.\n');
    exit(0);
end

% --- 辅助函数：根据名字获取函数句柄 ---
function func = get_drive_function(type_name)
    switch type_name
        case 'TARE'
            func = @TARE_spatial;
        case '2-locus_TARE'
            func = @Two_locus_TARE_spatial;
        case 'TADE_suppression'
            func = @TADE_suppression_spatial;
        case 'CifAB'
            func = @CifAB_spatial;
        case 'Wolbachia'
            func = @Wolbachia_spatial;
        otherwise
            error('Unknown drive type: %s', type_name);
    end
end

% --- 辅助函数：保存单个结果 ---
function save_single_result(result, drive_type, r)
    output_filename = sprintf('./results/variable_release/%s_optimal_variable_release.csv', drive_type);
    output_dir = fileparts(output_filename);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    header = {'time'};
    f_cols = arrayfun(@(k) sprintf('f%d', k), 1:r, 'UniformOutput', false);
    header = [header, f_cols, {'efficiency', 'exitflag'}];
    
    new_row = [result.T, result.x_final, -result.f_final, result.exitflag];
    new_table = array2table(new_row, 'VariableNames', header);
    
    if exist(output_filename, 'file')
        try
            existing = readtable(output_filename);
            existing(existing.time == result.T, :) = [];
            final = sortrows([existing; new_table], 'time');
            writetable(final, output_filename);
        catch
            writetable(new_table, [output_filename, '.new.csv']);
        end
    else
        writetable(new_table, output_filename);
    end
end