function [res, initial, final, rsize] = TARE_spatial(intro, Nr, T)
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

%% introduction
rr = length(intro);
A(4, 1:rr) = intro./(1-intro)*K; % dd introduction
initial = 2*A(4, :)*ring_area;
rsize = A(4, :)*ring_area;

%% reaction-diffusion
for t = 1:Nt

    %% time t
    N = sum(A, 1);
    fww = A(1, :)./ N;
    fwd = A(2, :)./ N;
    fwr = A(3, :)./ N;
    fdd = A(4, :)./ N;
    fdr = A(5, :)./ N;
    frr = A(6, :)./ N;

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
end
ndd = A(4, 1:end);
nwd = A(2, 1:end);
ndr = A(5, 1:end);
nd = 2*ndd+nwd+ndr;
final = nd*ring_area;
res = final / initial;
end