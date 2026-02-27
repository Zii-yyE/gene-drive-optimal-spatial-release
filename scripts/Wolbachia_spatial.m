function [res, initial, final, rsize] = Wolbachia_spatial(intro, Nr, T)
%% mesh generation
% Spatial grid
dr = 0.5;            % Spatial step size
ring_area = (1:2:(2*Nr-1))';
R = (Nr-1/2)*dr;
r = dr/2:dr:R;
inv_sqrdr = 1./(dr.^2);

% Temporal grid
dt = 0.05;          % Time step size 
Nt = T / dt;       % Number of time steps

%% initialization
lambda = 9;
K = 1;
% m: number of genotypes
m = 2;
A = zeros(m, Nr); 
A(1, :) = K;        % w population

%% introduction
rr = length(intro);
A(2, 1:rr) = intro./(1-intro)*K;
initial = A(2, :)*ring_area;
rsize = A(2, :)*ring_area;

%% reaction-diffusion
for t = 1:Nt

    %% time t
    N = sum(A, 1);
    fw = A(1, :)./ N;
    fd = A(2, :)./ N;
    
    %% reaction
    P = zeros(m, Nr);
    P(1, :) = fw.^2;
    P(2, :) = fd .* 0.75 .* (fw + fd .* 0.75);
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
nd = A(2, 1:end);
nw = A(1, 1:end);
final = nd*ring_area;
res = final / initial;
end