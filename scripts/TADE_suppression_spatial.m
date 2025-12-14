function [res, initial, final] = TADE_suppression_spatial(intro, Nr, T)
%% mesh generation
% Spatial grid
dr = 0.5;           
ring_area = (1:2:(2*Nr-1))';
R = (Nr-1/2)*dr;
r = dr/2:dr:R;
inv_sqrdr = 1./(dr.^2);

% Temporal grid
dt = 0.05;          
Nt = T / dt;      

%% initialization
lambda = 9;
K = 1;
m = 6;              % m: number of genotypes 
A = zeros(m, Nr); 
A(1, :) = K;        % ww population

%% beginning population size
pop_start = sum(A, 1)*ring_area;

%% introduction
rr = length(intro);
A(3, 1:rr) = intro./(1-intro)*K; % wdwr introduction
initial = A(3, :)*ring_area;

%% reaction-diffusion
for t = 1:Nt
    
    % time t
    N = sum(A, 1);
    fww_ww = A(1, :)./ N;
    fwd_ww = A(2, :)./ N;
    fwd_wr = A(3, :)./ N;
    fdd_ww = A(4, :)./ N;
    fdd_wr = A(5, :)./ N;
    fdd_rr = A(6, :)./ N;
    
    %% reaction
    P = zeros(6, Nr);
    P(1, :) = fww_ww.^2;
    P(3, :) = fww_ww.*(fwd_wr+fwd_ww+fdd_ww+fdd_wr+fdd_rr);
    P(6, :) = fwd_wr.*(0.25*fwd_wr+0.5*fwd_ww+0.5*fdd_rr+0.25*fdd_wr+0.5*fdd_ww) ...
        + fwd_ww.*(0.25*fwd_ww+0.5*fdd_rr+0.5*fdd_wr+0.5*fdd_ww);
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
end
pop_end = sum(A, 1)*ring_area;
pop_delta = pop_end - pop_start;
final = -pop_delta;
res = final / initial;
end