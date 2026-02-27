function run_all_circle_optimization()
    % All-in-One Circle Optimization Script with START POINTS
    % Uses Local Hill Climbing to drastically speed up optimization.
    % T=10000
    
    try
        % ================== 1. 任务配置 (请在此处输入你的 Start Points) ==================
        % 格式: {DriveType, StartRadius(整数), StartFrequency(0-1)}
        % 程序会从这个点开始局部搜索，不用跑完所有半径。
        jobs_config = {
            'TARE',             5, 0.09; 
            '2-locus_TARE',     10, 0.37;
            'TADE_suppression', 14, 0.06;
            'CifAB',            30, 0.45;
            'Wolbachia',        23, 0.41;
        };
        
        target_T = 10000; 
        R = 10000;        
        r_limit = 50;     % 半径搜索上限
        f_window = 0.15;  % 频率搜索窗口 (StartFreq +/- 0.15)
        
        fprintf('--- Starting Accelerated Circle Optimization (Local Search) ---\n');

        % ================== 2. 并行池设置 ==================
        delete(gcp('nocreate')); 
        parpool(length(jobs_config)); 
        
        fprintf('\n--- Parallel pool started. Launching workers... ---\n');
        
        % ================== 3. 并行执行所有任务 ==================
        parfor i = 1:size(jobs_config, 1)
            maxNumCompThreads(1);
            
            drive_type = jobs_config{i, 1};
            start_r    = jobs_config{i, 2};
            start_f    = jobs_config{i, 3};
            
            fprintf('[Worker %d] Optimizing %s (Start: r=%d, f=%.3f)...\n', i, drive_type, start_r, start_f);
            
            sim_func = get_drive_function(drive_type);
            
            % --- 局部爬山算法 (Local Hill Climbing for Radius) ---
            % 1. 评估起点
            current_r = max(1, min(r_limit, round(start_r)));
            [best_eff, best_f] = optimize_frequency_at_r(current_r, start_f, f_window, R, target_T, sim_func);
            
            while true
                % 尝试左右邻居
                neighbors = [current_r - 1, current_r + 1];
                improved = false;
                
                for next_r = neighbors
                    if next_r < 1 || next_r > r_limit, continue; end
                    
                    % 优化邻居的频率 (以当前最佳 f 为中心缩小范围)
                    [eff, f_opt] = optimize_frequency_at_r(next_r, best_f, f_window, R, target_T, sim_func);
                    
                    if eff > best_eff
                        best_eff = eff;
                        best_f = f_opt;
                        current_r = next_r;
                        improved = true;
                        break; % 贪婪策略：一发现更好的就移动
                    end
                end
                
                if ~improved
                    break; % 左右都不如中间，到达局部峰值，停止
                end
            end
            
            % --- 4. 结果保存 ---
            res_struct = struct('time', target_T, ...
                                'radius', current_r, ...
                                'frequency', best_f, ...
                                'efficiency', best_eff);
            
            save_single_result(res_struct, drive_type);
            
            fprintf('[Worker %d] Finished %s. Best: R=%d, F=%.3f, Eff=%.3f\n', ...
                    i, drive_type, current_r, best_f, best_eff);
        end
        
    catch ME
        fprintf('Global Error occurred:\n%s\n', ME.message);
        disp(ME.stack);
        exit(1);
    end
    
    fprintf('All jobs complete.\n');
    exit(0);
end

% --- 辅助函数：在固定半径下优化频率 ---
function [efficiency, optimal_f] = optimize_frequency_at_r(r, center_f, window, R, T, sim_func)
    % 缩小搜索范围: [center - window, center + window]
    lb = max(0, center_f - window);
    ub = min(1, center_f + window);
    
    obj_fun = @(f) -drive_efficiency_circle_helper(r, f, R, T, sim_func);
    
    % 精度 1e-3
    options = optimset('TolX', 1e-3, 'Display', 'off');
    
    [optimal_f, neg_eff] = fminbnd(obj_fun, lb, ub, options);
    efficiency = -neg_eff;
end

% --- 辅助函数：Circle 仿真包装器 ---
function output = drive_efficiency_circle_helper(radius, frequency, R, T, drive_function)
    intro_vec = zeros(1, R);
    r_idx = min(round(radius), R);
    if r_idx > 0
        intro_vec(1:r_idx) = frequency;
    end
    output = drive_function(intro_vec, R, T);
end

% --- 辅助函数：根据名字获取函数句柄 ---
function func = get_drive_function(type_name)
    switch type_name
        case 'TARE', func = @TARE_spatial;
        case '2-locus_TARE', func = @Two_locus_TARE_spatial;
        case 'TADE_suppression', func = @TADE_suppression_spatial;
        case 'CifAB', func = @CifAB_spatial;
        case 'Wolbachia', func = @Wolbachia_spatial;
        otherwise, error('Unknown drive type: %s', type_name);
    end
end

% --- 辅助函数：保存结果 ---
function save_single_result(result, drive_type)
    output_dir = './results/circle_release';
    output_filename = sprintf('%s/%s_optimal_circle_release.csv', output_dir, drive_type);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    new_table = struct2table(result);
    if exist(output_filename, 'file')
        try
            existing = readtable(output_filename);
            existing(existing.time == result.time, :) = [];
            final = sortrows([existing; new_table], 'time');
            writetable(final, output_filename);
        catch
            writetable(new_table, [output_filename, '.new.csv']);
        end
    else
        writetable(new_table, output_filename);
    end
end