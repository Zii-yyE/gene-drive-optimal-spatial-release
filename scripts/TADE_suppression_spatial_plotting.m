function res = TADE_suppression_spatial_plotting(intro, Nr, T, varargin)
% TADE suppression model with real-time plotting of population size.
% intro: initial frequency of drive individuals
% Nr:    Domain radius
% T:     Simulation time.
% plot_update_frequency (optional): Update plot every Nth dt step.

%% Plotting Setup
persistent plot_fig_handle plot_lines_handles r_plot_persistent

enable_plotting_in_loop = true;
if nargin > 3 && ~isempty(varargin{1})
    plot_update_dts = varargin{1};
    if plot_update_dts <= 0 || isinf(plot_update_dts)
        enable_plotting_in_loop = false;
    end
else
    plot_update_dts = 20;
end

if enable_plotting_in_loop
    if isempty(r_plot_persistent) || length(r_plot_persistent) ~= Nr
        plot_fig_handle = [];
        plot_lines_handles = [];
    end
end

%% Mesh Generation
dr = 0.5;
ring_area = (1:2:(2*Nr-1))';
R = (Nr-1/2)*dr;
r = dr/2:dr:R;
inv_sqrdr = 1./(dr.^2);
inv_rdr = 1./(dr*r);
dt = 0.05;
Nt = T / dt;

%% Initialization
lambda = 9;
K = 1;
m = 6;
A = zeros(m, Nr);
A(1, :) = K; % ww population
genotype_names = {'ww-ww', 'wd-ww', 'wd-wr', 'dd-ww', 'dd-wr', 'dd-rr'};

%% Initial Population Size
pop_start = sum(A, 1)*ring_area;

%% Introduction
rr = length(intro);
A(3, 1:rr) = intro./(1-intro)*K; % wd-wr introduction
released = A(3, :)*ring_area;
fprintf('\n');
fprintf('Number of drive individuals released: %f\n', released);


%% Setup Comparison Plot (Initial vs. Final State)
figure('Name', 'Population Size: Initial vs. Final State');
ax_left = subplot(1, 2, 1);
ax_right = subplot(1, 2, 2);

hold(ax_left, 'on');
colors = lines(m);
for k_geno = 1:m
    plot(ax_left, r, A(k_geno,:), 'LineWidth', 1.5, 'Color', colors(k_geno,:), 'DisplayName', genotype_names{k_geno});
end
title(ax_left, 'Initial State (t=0)');
xlabel(ax_left, 'Distance to Origin');
ylabel(ax_left, 'Population Size');
grid(ax_left, 'on');
pbaspect(ax_left, [1 1 1]);
hold(ax_left, 'off');

title(ax_right, sprintf('Final State (t=%d)', T));
xlabel(ax_right, 'Distance to Origin');
grid(ax_right, 'on');
pbaspect(ax_right, [1 1 1]);
set(ax_right, 'YTickLabel', []);


%% Setup Real-time Plot
if enable_plotting_in_loop && (isempty(plot_fig_handle) || ~isvalid(plot_fig_handle))
    plot_fig_handle = figure('Name', 'Real-Time Population Size');
    ax_plot = axes(plot_fig_handle);
    hold(ax_plot, 'on');
    plot_lines_handles = gobjects(m,1);
    for k_geno = 1:m
        plot_lines_handles(k_geno) = plot(ax_plot, r, A(k_geno,:), 'LineWidth', 1.5, 'Color', colors(k_geno,:), 'DisplayName', genotype_names{k_geno});
    end
    xlabel(ax_plot, 'Distance to Origin');
    ylabel(ax_plot, 'Population Size');
    grid(ax_plot, 'on');
    legend(ax_plot, 'show', 'Location', 'eastoutside');
    r_plot_persistent = r;
end

%% Reaction-Diffusion Loop
for t = 1:Nt
    % Time t
    N = sum(A, 1);
    fww_ww = A(1, :)./ N;
    fwd_ww = A(2, :)./ N;
    fwd_wr = A(3, :)./ N;
    fdd_ww = A(4, :)./ N;
    fdd_wr = A(5, :)./ N;
    fdd_rr = A(6, :)./ N;

    % Reaction
    P = zeros(6, Nr);
    P(1, :) = fww_ww.^2;
    P(3, :) = fww_ww.*(fwd_wr+fwd_ww+fdd_ww+fdd_wr+fdd_rr);
    P(6, :) = fwd_wr.*(0.25*fwd_wr+0.5*fwd_ww+0.5*fdd_rr+0.25*fdd_wr+0.5*fdd_ww) ...
        + fwd_ww.*(0.25*fwd_ww+0.5*fdd_rr+0.5*fdd_wr+0.5*fdd_ww);
    rxn = lambda./((lambda-1)*N+1);
    reaction = rxn.*N.*P-A;

    % Diffusion
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

    % Update
    A = A + dt * (reaction + diffusion);
    A = max(A, 0);

    % Real-time Plotting
    if enable_plotting_in_loop && (mod(t - 1, plot_update_dts) == 0 || t == Nt)
        if isvalid(plot_fig_handle)
            for k_geno = 1:m
                set(plot_lines_handles(k_geno), 'YData', A(k_geno,:));
            end
            title(plot_fig_handle.CurrentAxes, sprintf('Population Size at Time = %.2f', t*dt));
            drawnow;
            pause(0.0001);
        end
    end
end

%% Final Population and Suppression Calculation
pop_end = sum(A, 1)*ring_area;
pop_delta = pop_end - pop_start;
suppressed = -pop_delta;
fprintf('Number of individuals suppressed: %f\n', suppressed);
res = suppressed / released;

%% Plot Final State on Comparison Figure
hold(ax_right, 'on');
for k_geno = 1:m
    plot(ax_right, r, A(k_geno,:), 'LineWidth', 1.5, 'Color', colors(k_geno,:), 'DisplayName', genotype_names{k_geno});
end
legend(ax_right, 'show', 'Location', 'eastoutside');
hold(ax_right, 'off');

end