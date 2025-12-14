function [res, initial, final] = CifAB_spatial(intro, Nr, T)
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
m = 3;
A = zeros(m, Nr); 
A(1, :) = K;        % ww population

%% introduction
rr = length(intro);
A(3, 1:rr) = intro./(1-intro)*K;
initial = 2*A(3, :)*ring_area;

%% reaction-diffusion
for t = 1:Nt

    %% time t
    N = sum(A, 1);
    fww = A(1, :)./ N;
    fwd = A(2, :)./ N;
    fdd = A(3, :)./ N;
    
    %% reaction
    P = zeros(3, Nr);
    P(1, :) = fww.^2+0.5.*fww.*fwd+0.25.*fwd.^2;
    P(2, :) = fww.*fdd+0.5.*fww.*fwd+0.5.*fwd.^2+fwd.*fdd;
    P(3, :) = fdd.^2+fdd.*fwd+0.25.*fwd.^2;
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
ndd = A(3, 1:end);
nwd = A(2, 1:end);
final = (2*ndd+nwd)*ring_area;
res = final / initial;
end