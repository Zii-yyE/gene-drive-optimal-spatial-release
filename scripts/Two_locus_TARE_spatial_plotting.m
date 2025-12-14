function res = Two_Locus_TARE_spatial_plotting(intro, Nr, T, varargin)
%   Two-locus TARE reaction-diffusion model with real-time plotting.
%   intro: initial frequency of drive individuals (genotype dddd)
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
m = 25;             % m: number of genotypes
A = zeros(m, Nr);
A(1, :) = K;        % ww-ww population

% Generate genotype names programmatically for the 5x5 system
g_alleles = {'ww', 'wd', 'wr', 'dd', 'dr'};
genotype_names = cell(m, 1);
idx = 1;
for i = 1:5
    for j = 1:5
        genotype_names{idx} = sprintf('%s-%s', g_alleles{i}, g_alleles{j});
        idx = idx + 1;
    end
end
allele_names = {'w allele', 'd allele', 'r allele'};

%% Introduction
rr = length(intro);
A(19, 1:rr) = intro./(1-intro)*K; % dd-dd introduction
released = 2*A(19, :)*ring_area;
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
colors_g = jet(m); % FIX: Use 'jet' colormap for 25 unique colors
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

% Calculate initial allele frequencies (using first locus)
F_ww = sum(F_initial_g(1:5,:),1); F_wd = sum(F_initial_g(6:10,:),1); F_wr = sum(F_initial_g(11:15,:),1);
F_dd = sum(F_initial_g(16:20,:),1); F_dr = sum(F_initial_g(21:25,:),1);
F_initial_a = [F_ww + 0.5*(F_wd + F_wr); F_dd + 0.5*(F_wd + F_dr); 0.5*(F_wr + F_dr)];

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
    fwwww = A(1, :)./N; fwwwd = A(2, :)./N; fwwwr = A(3, :)./N; fwwdd = A(4, :)./N; fwwdr = A(5, :)./N;
    fwdww = A(6, :)./N; fwdwd = A(7, :)./N; fwdwr = A(8, :)./N; fwddd = A(9, :)./N; fwddr = A(10, :)./N;
    fwrww = A(11, :)./N; fwrwd = A(12, :)./N; fwrwr = A(13, :)./N; fwrdd = A(14, :)./N; fwrdr = A(15, :)./N;
    fddww = A(16, :)./N; fddwd = A(17, :)./N; fddwr = A(18, :)./N; fdddd = A(19, :)./N; fdddr = A(20, :)./N;
    fdrww = A(21, :)./N; fdrwd = A(22, :)./N; fdrwr = A(23, :)./N; fdrdd = A(24, :)./N; fdrdr = A(25, :)./N;

    %% reaction
    P = zeros(25, Nr);
    P(1, :) = fwwww.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwwwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(2, :) = 0;
    P(3, :) = fwwww.*(0.5*fwwwr+0.25*fwrwr)+fwwwr.*(0.5*fwwww+0.5*fwwwr+0.25*fwrww+0.25*fwrwr)+fwdww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrww.*(0.25*fwwwr+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.25*fwwwr+0.125*fwrww+0.125*fwrwr);
    P(4, :) = 0; P(5, :) = 0; P(6, :) = 0;
    P(7, :) = fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwddd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fddwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdddd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fdddr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrdd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(8, :) = fwdww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwdwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fddww.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fddwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fddwr.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fdddr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(9, :) = 0;
    P(10, :) = fwdwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwddd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwddr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fddwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdddd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fdddr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fdrdd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrdr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr);
    P(11, :) = fwwww.*(0.5*fwrww+0.25*fwrwr)+fwwwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwwwr.*(0.25*fwrww+0.125*fwrwr)+fwrww.*(0.5*fwwww+0.25*fwwwr+0.5*fwrww+0.25*fwrwr)+fwrwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.125*fwwwr+0.25*fwrww+0.125*fwrwr);
    P(12, :) = fwwwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwwdd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwwdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwddd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwrwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrdd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwrdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrdd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(13, :) = fwwww.*(0.25*fwrwr)+fwwwd.*(0.25*fwwwr+0.125*fwrwr)+fwwwr.*(0.25*fwrww+0.25*fwrwr)+fwwdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdww.*(0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwdwr.*(0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwrww.*(0.25*fwwwr+0.25*fwrwr)+fwrwd.*(0.25*fwwwr+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.25*fwwwr+0.25*fwrww+0.25*fwrwr)+fwrdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(14, :) = 0;
    P(15, :) = fwwwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwwdd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fwwdr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwdwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwddd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwddr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwrwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwrdd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fwrdr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fdrdd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrdr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr);
    P(16, :) = 0; P(17, :) = 0; P(18, :) = 0;
    P(19, :) = fwdwd.*(0.0625*fwdwd+0.125*fwddd+0.0625*fwddr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.0625*fdrwd+0.125*fdrdd+0.0625*fdrdr)+fwddd.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwddr.*(0.0625*fwdwd+0.125*fwddd+0.0625*fwddr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.0625*fdrwd+0.125*fdrdd+0.0625*fdrdr)+fddwd.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdddd.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdddr.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrwd.*(0.0625*fwdwd+0.125*fwddd+0.0625*fwddr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.0625*fdrwd+0.125*fdrdd+0.0625*fdrdr)+fdrdd.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrdr.*(0.0625*fwdwd+0.125*fwddd+0.0625*fwddr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.0625*fdrwd+0.125*fdrdd+0.0625*fdrdr);
    P(20, :) = fwdww.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwdwd.*(0.125*fwdww+0.125*fwdwd+0.125*fwdwr+0.125*fwddd+0.125*fwddr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.125*fdrww+0.125*fdrwd+0.125*fdrwr+0.125*fdrdd+0.125*fdrdr)+fwdwr.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwddd.*(0.25*fwdww+0.125*fwdwd+0.25*fwdwr+0.125*fwddr+0.5*fddww+0.25*fddwd+0.5*fddwr+0.25*fdddr+0.25*fdrww+0.125*fdrwd+0.25*fdrwr+0.125*fdrdr)+fwddr.*(0.125*fwdww+0.125*fwdwd+0.125*fwdwr+0.125*fwddd+0.125*fwddr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.125*fdrww+0.125*fdrwd+0.125*fdrwr+0.125*fdrdd+0.125*fdrdr)+fddww.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fddwd.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fddwr.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdddd.*(0.5*fwdww+0.25*fwdwd+0.5*fwdwr+0.25*fwddr+1.0*fddww+0.5*fddwd+1.0*fddwr+0.5*fdddr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fdddr.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fdrww.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrwd.*(0.125*fwdww+0.125*fwdwd+0.125*fwdwr+0.125*fwddd+0.125*fwddr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.125*fdrww+0.125*fdrwd+0.125*fdrwr+0.125*fdrdd+0.125*fdrdr)+fdrwr.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrdd.*(0.25*fwdww+0.125*fwdwd+0.25*fwdwr+0.125*fwddr+0.5*fddww+0.25*fddwd+0.5*fddwr+0.25*fdddr+0.25*fdrww+0.125*fdrwd+0.25*fdrwr+0.125*fdrdr)+fdrdr.*(0.125*fwdww+0.125*fwdwd+0.125*fwdwr+0.125*fwddd+0.125*fwddr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.125*fdrww+0.125*fdrwd+0.125*fdrwr+0.125*fdrdd+0.125*fdrdr);
    P(21, :) = 0;
    P(22, :) = fwdwd.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fwddd.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fwddr.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fddwd.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdddd.*(0.5*fwwwd+0.5*fwrww+0.5*fwrwd+0.25*fwrwr)+fdddr.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdrwd.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fdrdd.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdrdr.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr);
    P(23, :) = fwdww.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fwdwd.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fwdwr.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fwddr.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fddww.*(0.5*fwwwd+0.5*fwrww+0.5*fwrwd+0.25*fwrwr)+fddwd.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fddwr.*(0.5*fwwwd+0.5*fwrww+0.5*fwrwd+0.25*fwrwr)+fdddr.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdrww.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdrwd.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr)+fdrwr.*(0.25*fwwwd+0.25*fwrww+0.25*fwrwd+0.125*fwrwr)+fdrdr.*(0.125*fwwwd+0.125*fwrww+0.125*fwrwd+0.0625*fwrwr);
    P(24, :) = fwwwd.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwwdd.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwwdr.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwdwd.*(0.125*fwwwd+0.25*fwwdd+0.125*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.125*fwrwd+0.25*fwrdd+0.125*fwrdr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwddd.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwddr.*(0.125*fwwwd+0.25*fwwdd+0.125*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.125*fwrwd+0.25*fwrdd+0.125*fwrdr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwrwd.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fwrdd.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwrdr.*(0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fddwd.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdddd.*(0.5*fwwwd+1.0*fwwdd+0.5*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fwrwd+1.0*fwrdd+0.5*fwrdr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdddr.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrwd.*(0.125*fwwwd+0.25*fwwdd+0.125*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.125*fwrwd+0.25*fwrdd+0.125*fwrdr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr)+fdrdd.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdrdr.*(0.125*fwwwd+0.25*fwwdd+0.125*fwwdr+0.125*fwdwd+0.25*fwddd+0.125*fwddr+0.125*fwrwd+0.25*fwrdd+0.125*fwrdr+0.125*fddwd+0.25*fdddd+0.125*fdddr+0.125*fdrwd+0.25*fdrdd+0.125*fdrdr);
    P(25, :) = fwwww.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwwwd.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fwwwr.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwwdd.*(0.5*fwdww+0.25*fwdwd+0.5*fwdwr+0.25*fwddr+1.0*fddww+0.5*fddwd+1.0*fddwr+0.5*fdddr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fwwdr.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fwdww.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwdwd.*(0.125*fwwwd+0.25*fwwdd+0.25*fwwdr+0.125*fwdww+0.25*fwdwd+0.125*fwdwr+0.25*fwddd+0.25*fwddr+0.125*fwrwd+0.0625*fwrwr+0.25*fwrdd+0.25*fwrdr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fwdwr.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwddd.*(0.25*fwwdr+0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddr+0.125*fwrwr+0.25*fwrdr+0.5*fddww+0.25*fddwd+0.5*fddwr+0.25*fdddr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fwddr.*(0.125*fwwwd+0.25*fwwdd+0.25*fwwdr+0.125*fwdww+0.25*fwdwd+0.125*fwdwr+0.25*fwddd+0.25*fwddr+0.125*fwrwd+0.0625*fwrwr+0.25*fwrdd+0.25*fwrdr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fwrww.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwrwd.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fwrwr.*(0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fddwd+1.0*fdddd+0.5*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fwrdd.*(0.5*fwdww+0.25*fwdwd+0.5*fwdwr+0.25*fwddr+1.0*fddww+0.5*fddwd+1.0*fddwr+0.5*fdddr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fwrdr.*(0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddd+0.25*fwddr+0.5*fddww+0.5*fddwd+0.5*fddwr+0.5*fdddd+0.5*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fddww.*(0.5*fwwwd+1.0*fwwdd+0.5*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fwrwd+1.0*fwrdd+0.5*fwrdr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fddwd.*(0.25*fwwwd+0.5*fwwdd+0.5*fwwdr+0.25*fwdwd+0.25*fwddd+0.25*fwddr+0.25*fwrwd+0.125*fwrwr+0.5*fwrdd+0.5*fwrdr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fddwr.*(0.5*fwwwd+1.0*fwwdd+0.5*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.5*fwrwd+1.0*fwrdd+0.5*fwrdr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdddd.*(0.5*fwwdr+0.25*fwdwd+0.25*fwddr+0.25*fwrwr+0.5*fwrdr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fdddr.*(0.25*fwwwd+0.5*fwwdd+0.5*fwwdr+0.25*fwdwd+0.25*fwddd+0.25*fwddr+0.25*fwrwd+0.125*fwrwr+0.5*fwrdd+0.5*fwrdr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fdrww.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdrwd.*(0.125*fwwwd+0.25*fwwdd+0.25*fwwdr+0.125*fwdww+0.25*fwdwd+0.125*fwdwr+0.25*fwddd+0.25*fwddr+0.125*fwrwd+0.0625*fwrwr+0.25*fwrdd+0.25*fwrdr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr)+fdrwr.*(0.25*fwwwd+0.5*fwwdd+0.25*fwwdr+0.25*fwdwd+0.5*fwddd+0.25*fwddr+0.25*fwrwd+0.5*fwrdd+0.25*fwrdr+0.25*fddwd+0.5*fdddd+0.25*fdddr+0.25*fdrwd+0.5*fdrdd+0.25*fdrdr)+fdrdd.*(0.25*fwwdr+0.25*fwdww+0.25*fwdwd+0.25*fwdwr+0.25*fwddr+0.125*fwrwr+0.25*fwrdr+0.5*fddww+0.25*fddwd+0.5*fddwr+0.25*fdddr+0.5*fdrww+0.25*fdrwd+0.5*fdrwr+0.25*fdrdr)+fdrdr.*(0.125*fwwwd+0.25*fwwdd+0.25*fwwdr+0.125*fwdww+0.25*fwdwd+0.125*fwdwr+0.25*fwddd+0.25*fwddr+0.125*fwrwd+0.0625*fwrwr+0.25*fwrdd+0.25*fwrdr+0.25*fddww+0.25*fddwd+0.25*fddwr+0.25*fdddd+0.25*fdddr+0.25*fdrww+0.25*fdrwd+0.25*fdrwr+0.25*fdrdd+0.25*fdrdr);
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

    %% Real-time Plotting (inside the time loop)
    if enable_plotting_in_loop && (mod(t - 1, plot_update_dts) == 0 || t == Nt)
        N_plot = sum(A,1);
        N_plot_eff = max(N_plot, 1e-9);
        F_plot_g = A ./ N_plot_eff;

        % Calculate allele frequencies for plotting
        F_ww = sum(F_plot_g(1:5,:),1); F_wd = sum(F_plot_g(6:10,:),1); F_wr = sum(F_plot_g(11:15,:),1);
        F_dd = sum(F_plot_g(16:20,:),1); F_dr = sum(F_plot_g(21:25,:),1);
        F_plot_a = [F_ww + 0.5*(F_wd + F_wr); F_dd + 0.5*(F_wd + F_dr); 0.5*(F_wr + F_dr)];

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
        drawnow limitrate;
    end
end

%% Plot Final States on Comparison Figures
% Final Genotype Frequencies
N_final = sum(A, 1);
F_final_g = A ./ max(N_final, 1e-9);
hold(ax_right_g, 'on'); % FIX: Explicitly set hold on for the final plot
for k_geno = 1:m
    plot(ax_right_g, r, F_final_g(k_geno,:), 'LineWidth', 1.5, 'Color', colors_g(k_geno,:), 'DisplayName', genotype_names{k_geno});
end
legend(ax_right_g, 'show', 'Location', 'eastoutside');
hold(ax_right_g, 'off');

% Final Allele Frequencies
F_ww = sum(F_final_g(1:5,:),1); F_wd = sum(F_final_g(6:10,:),1); F_wr = sum(F_final_g(11:15,:),1);
F_dd = sum(F_final_g(16:20,:),1); F_dr = sum(F_final_g(21:25,:),1);
F_final_a = [F_ww + 0.5*(F_wd + F_wr); F_dd + 0.5*(F_wd + F_dr); 0.5*(F_wr + F_dr)];
hold(ax_right_a, 'on'); % FIX: Explicitly set hold on for the final plot
for k_allele = 1:3
    plot(ax_right_a, r, F_final_a(k_allele,:), 'LineWidth', 2, 'Color', colors_a(k_allele,:), 'DisplayName', allele_names{k_allele});
end
legend(ax_right_a, 'show', 'Location', 'northeast');
hold(ax_right_a, 'off');

%% Calculate Harvested Individuals
nd = 2*sum(A(16:20,:),1) + sum(A(6:10,:),1) + sum(A(21:25,:),1);
harvested = nd*ring_area;
fprintf('Number of drive alleles harvested: %f\n', harvested);
res = harvested / released;

end