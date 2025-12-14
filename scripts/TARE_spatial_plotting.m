function res = TARE_spatial_plotting(intro, Nr, T, varargin)
%   TARE reaction-diffusion model with real-time plotting for genotypes and alleles.
%   intro: initial frequency of drive individuals (genotype dd)
%   Nr:    Domain radius
%   T:     Simulation time.
%   plot_update_frequency (optional): Update plot every Nth dt step. Default is 20.
%                                     Set to 0 or Inf to disable plotting within the loop.

%% Plotting Setup
% Persistent handles for the real-time animation windows
persistent plot_fig_handle_g plot_lines_handles_g r_plot_persistent_g
persistent plot_fig_handle_a plot_lines_handles_a r_plot_persistent_a

enable_plotting_in_loop = true;
if nargin > 3 && ~isempty(varargin{1})
    plot_update_dts = varargin{1};
    if plot_update_dts <= 0 || isinf(plot_update_dts)
        enable_plotting_in_loop = false;
    end
else
    plot_update_dts = 20; % Default: update plot every 20 dt steps
end

if enable_plotting_in_loop
    % Clear persistent handles if Nr changes or for a fresh plot start
    if isempty(r_plot_persistent_g) || length(r_plot_persistent_g) ~= Nr
        plot_fig_handle_g = []; plot_lines_handles_g = [];
        plot_fig_handle_a = []; plot_lines_handles_a = [];
    end
end

%% Mesh Generation
% Spatial grid
dr = 0.5;
ring_area = (1:2:(2*Nr-1))';
R = (Nr-1/2)*dr;
r = dr/2:dr:R;
inv_sqrdr = 1./(dr.^2);
inv_rdr = 1./(dr*r);
% Temporal grid
dt = 0.05;
Nt = T / dt;

%% Initialization
lambda = 9;
K = 1;
m = 6;
A = zeros(m, Nr);
A(1, :) = K;        % ww population
genotype_names = {'ww', 'wd', 'wr', 'dd', 'dr', 'rr'};
allele_names = {'w allele', 'd allele', 'r allele'};

%% Introduction
rr = length(intro);
A(4, 1:rr) = intro./(1-intro)*K; % dd introduction
released = 2*A(4, :)*ring_area;
fprintf('\n');
fprintf('Number of drive alleles released: %f\n', released);

%% Setup Comparison Plots (Initial and Final States)

% --- Figure 1: Genotype Frequencies ---
figure('Name', 'Genotype Frequencies: Initial vs. Final State');
ax_left_g = subplot(1, 2, 1);
ax_right_g = subplot(1, 2, 2);
N_initial = sum(A, 1);
F_initial_g = A ./ max(N_initial, 1e-9);

hold(ax_left_g, 'on');
colors_g = lines(m);
for k_geno = 1:m
    plot(ax_left_g, r, F_initial_g(k_geno,:), 'LineWidth', 1.5, 'Color', colors_g(k_geno,:), 'DisplayName', genotype_names{k_geno});
end
title(ax_left_g, 'Initial State (t=0)');
xlabel(ax_left_g, 'Distance to Origin');
ylabel(ax_left_g, 'Genotype Frequency');
grid(ax_left_g, 'on'); ylim(ax_left_g, [-0.05, 1.05]); pbaspect(ax_left_g, [1 1 1]);
hold(ax_left_g, 'off');

title(ax_right_g, sprintf('Final State (t=%d)', T));
xlabel(ax_right_g, 'Distance to Origin');
grid(ax_right_g, 'on'); ylim(ax_right_g, [-0.05, 1.05]); pbaspect(ax_right_g, [1 1 1]);
set(ax_right_g, 'YTickLabel', []);

% --- Figure 2: Allele Frequencies ---
figure('Name', 'Allele Frequencies: Initial vs. Final State');
ax_left_a = subplot(1, 2, 1);
ax_right_a = subplot(1, 2, 2);

% Calculate initial allele frequencies
F_initial_a = [F_initial_g(1,:) + 0.5*(F_initial_g(2,:) + F_initial_g(3,:)); ... % w
               F_initial_g(4,:) + 0.5*(F_initial_g(2,:) + F_initial_g(5,:)); ... % d
               F_initial_g(6,:) + 0.5*(F_initial_g(3,:) + F_initial_g(5,:))];    % r

hold(ax_left_a, 'on');
colors_a = [0 0 1; 1 0 0; 0 1 0]; % Blue, Red, Green for w, d, r
for k_allele = 1:3
    plot(ax_left_a, r, F_initial_a(k_allele,:), 'LineWidth', 2, 'Color', colors_a(k_allele,:), 'DisplayName', allele_names{k_allele});
end
title(ax_left_a, 'Initial State (t=0)');
xlabel(ax_left_a, 'Distance to Origin');
ylabel(ax_left_a, 'Allele Frequency');
grid(ax_left_a, 'on'); ylim(ax_left_a, [-0.05, 1.05]); pbaspect(ax_left_a, [1 1 1]);
hold(ax_left_a, 'off');

title(ax_right_a, sprintf('Final State (t=%d)', T));
xlabel(ax_right_a, 'Distance to Origin');
grid(ax_right_a, 'on'); ylim(ax_right_a, [-0.05, 1.05]); pbaspect(ax_right_a, [1 1 1]);
set(ax_right_a, 'YTickLabel', []);

%% Setup Real-time Plots
if enable_plotting_in_loop
    % Window for Genotype real-time plot
    if isempty(plot_fig_handle_g) || ~isvalid(plot_fig_handle_g)
        plot_fig_handle_g = figure('Name', 'Real-Time Genotype Frequencies');
        ax_plot_g = axes(plot_fig_handle_g);
        hold(ax_plot_g, 'on');
        plot_lines_handles_g = gobjects(m,1);
        for k_geno = 1:m
            plot_lines_handles_g(k_geno) = plot(ax_plot_g, r, F_initial_g(k_geno,:), 'LineWidth', 1.5, 'Color', colors_g(k_geno,:), 'DisplayName', genotype_names{k_geno});
        end
        xlabel(ax_plot_g, 'Distance to Origin'); ylabel(ax_plot_g, 'Genotype Frequency');
        grid(ax_plot_g, 'on'); legend(ax_plot_g, 'show', 'Location', 'eastoutside');
        ylim(ax_plot_g, [-0.05, 1.05]);
        r_plot_persistent_g = r;
    end
    % Window for Allele real-time plot
    if isempty(plot_fig_handle_a) || ~isvalid(plot_fig_handle_a)
        plot_fig_handle_a = figure('Name', 'Real-Time Allele Frequencies');
        ax_plot_a = axes(plot_fig_handle_a);
        hold(ax_plot_a, 'on');
        plot_lines_handles_a = gobjects(3,1);
        for k_allele = 1:3
            plot_lines_handles_a(k_allele) = plot(ax_plot_a, r, F_initial_a(k_allele,:), 'LineWidth', 2, 'Color', colors_a(k_allele,:), 'DisplayName', allele_names{k_allele});
        end
        xlabel(ax_plot_a, 'Distance to Origin'); ylabel(ax_plot_a, 'Allele Frequency');
        grid(ax_plot_a, 'on'); legend(ax_plot_a, 'show', 'Location', 'northeast');
        ylim(ax_plot_a, [-0.05, 1.05]);
        r_plot_persistent_a = r;
    end
end

%% Reaction-Diffusion Loop
for t = 1:Nt
    %% time t
    N = sum(A, 1);
    fww = A(1, :)./ N; fwd = A(2, :)./ N; fwr = A(3, :)./ N;
    fdd = A(4, :)./ N; fdr = A(5, :)./ N; frr = A(6, :)./ N;

    %% reaction
    P = zeros(6, Nr);
    P(1, :) = fww.*(fww+0.5*fwr)+fwr.*(0.5*fww+0.25*fwr);
    P(2, :) = fwd.*(0.5*fww+0.25*fwr)+fdd.*(fww+0.5*fwr)+fdr.*(0.5*fww+0.25*fwr);
    P(3, :) = fww.*(0.5*fwr+frr)+fwd.*(0.5*fww+0.25*fwr)+fwr.*(0.5*fww+0.5*fwr+0.5*frr)+fdr.*(0.5*fww+0.25*fwr)+frr.*(fww+0.5*fwr);
    P(4, :) = fwd.*(0.25*fwd+0.5*fdd+0.25*fdr)+fdd.*(0.5*fwd+fdd+0.5*fdr)+fdr.*(0.25*fwd+0.5*fdd+0.25*fdr);
    P(5, :) = fww.*(0.5*fwd+1.0*fdd+0.5*fdr)+fwd.*(0.5*fwd+0.25*fwr+0.5*fdd+0.5*fdr+0.5*frr)+fwr.*(0.5*fwd+1.0*fdd+0.5*fdr)+fdd.*(0.5*fwd+0.5*fwr+0.5*fdr+1.0*frr)+fdr.*(0.5*fwd+0.25*fwr+0.5*fdd+0.5*fdr+0.5*frr)+frr.*(0.5*fwd+1.0*fdd+0.5*fdr);
    P(6, :) = fww.*(0.5*fwd+0.5*fdr)+fwd.*(0.25*fwd+0.25*fwr+0.25*fdr+0.5*frr)+fwr.*(0.5*fwd+0.25*fwr+0.5*fdr+0.5*frr)+fdr.*(0.25*fwd+0.25*fwr+0.25*fdr+0.5*frr)+frr.*(0.5*fwd+0.5*fwr+0.5*fdr+1.0*frr);
    rxn = lambda./((lambda-1)*N+1);
    reaction = rxn.*N.*P-A;

    %% diffusion
    laplacian = zeros(size(A));
    for i = 2:(Nr-1)
        ord2_diff = (A(:,i+1) - 2*A(:,i) + A(:,i-1)) * inv_sqrdr;
        ord1_diff = (A(:,i+1) - A(:,i-1)) / (2*dr);
        laplacian(:,i) = ord2_diff + (1/r(i)) * ord1_diff;
    end
    ord2_diff_origin = (A(:,2) - A(:,1)) * inv_sqrdr;
    ord1_diff_origin = (A(:,2) - A(:,1)) / (2*dr);
    laplacian(:,1) = ord2_diff_origin + (1/r(1)) * ord1_diff_origin;
    ord2_diff_outer = (A(:,Nr-1) - A(:,Nr)) * inv_sqrdr;
    ord1_diff_outer = (A(:,Nr) - A(:,Nr-1)) / (2*dr);
    laplacian(:,end) = ord2_diff_outer + (1/r(end)) * ord1_diff_outer;
    diffusion = laplacian;

    %% update
    A = A + dt * (reaction + diffusion);
    A = max(A, 0);
    A(6, :) = 0;

    %% Real-time Plotting (inside the time loop)
    if enable_plotting_in_loop && (mod(t - 1, plot_update_dts) == 0 || t == Nt)
        N_plot = sum(A,1);
        N_plot_eff = max(N_plot, 1e-9);
        F_plot_g = A ./ N_plot_eff;

        % Calculate allele frequencies for plotting
        F_plot_a = [F_plot_g(1,:) + 0.5*(F_plot_g(2,:) + F_plot_g(3,:)); ... % w
                    F_plot_g(4,:) + 0.5*(F_plot_g(2,:) + F_plot_g(5,:)); ... % d
                    F_plot_g(6,:) + 0.5*(F_plot_g(3,:) + F_plot_g(5,:))];    % r

        % Update Genotype Plot
        if isvalid(plot_fig_handle_g)
            for k_geno = 1:m
                set(plot_lines_handles_g(k_geno), 'YData', F_plot_g(k_geno,:));
            end
            title(plot_fig_handle_g.CurrentAxes, sprintf('Genotype Frequencies at Time = %.2f', t*dt));
        end

        % Update Allele Plot
        if isvalid(plot_fig_handle_a)
            for k_allele = 1:3
                set(plot_lines_handles_a(k_allele), 'YData', F_plot_a(k_allele,:));
            end
            title(plot_fig_handle_a.CurrentAxes, sprintf('Allele Frequencies at Time = %.2f', t*dt));
        end
        drawnow;
    end
end

%% Plot Final States on Comparison Figures
% Final Genotype Frequencies
N_final = sum(A, 1);
F_final_g = A ./ max(N_final, 1e-9);
hold(ax_right_g, 'on');
for k_geno = 1:m
    plot(ax_right_g, r, F_final_g(k_geno,:), 'LineWidth', 1.5, 'Color', colors_g(k_geno,:), 'DisplayName', genotype_names{k_geno});
end
legend(ax_right_g, 'show', 'Location', 'eastoutside');
hold(ax_right_g, 'off');

% Final Allele Frequencies
F_final_a = [F_final_g(1,:) + 0.5*(F_final_g(2,:) + F_final_g(3,:)); ... % w
             F_final_g(4,:) + 0.5*(F_final_g(2,:) + F_final_g(5,:)); ... % d
             F_final_g(6,:) + 0.5*(F_final_g(3,:) + F_final_g(5,:))];    % r
hold(ax_right_a, 'on');
for k_allele = 1:3
    plot(ax_right_a, r, F_final_a(k_allele,:), 'LineWidth', 2, 'Color', colors_a(k_allele,:), 'DisplayName', allele_names{k_allele});
end
legend(ax_right_a, 'show', 'Location', 'northeast');
hold(ax_right_a, 'off');

%% Calculate Harvested Individuals
ndd = A(4, :);
nwd = A(2, :);
ndr = A(5, :);
harvested = (2*ndd+nwd+ndr)*ring_area;
fprintf('Number of drive alleles harvested: %f\n', harvested);
res = harvested / released;
end