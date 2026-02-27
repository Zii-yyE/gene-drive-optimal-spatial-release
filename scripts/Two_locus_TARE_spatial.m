function [res, initial, final, rsize] = Two_locus_TARE_spatial(intro, Nr, T)
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
m = 25;             % m: number of genotypes
A = zeros(m, Nr); 
A(1, :) = K;        % ww population

%% introduction
rr = length(intro);
A(19, 1:rr) = intro./(1-intro)*K; % dddd introduction
initial = 2*A(19, :)*ring_area;
rsize = A(19, :)*ring_area;

%% reaction-diffusion
for t = 1:Nt
    
    %% time t
    N = sum(A, 1);
    fwwww = A(1, :)./N;
    fwwwd = A(2, :)./N;
    fwwwr = A(3, :)./N;
    fwwdd = A(4, :)./N;
    fwwdr = A(5, :)./N;
    fwdww = A(6, :)./N;
    fwdwd = A(7, :)./N;
    fwdwr = A(8, :)./N;
    fwddd = A(9, :)./N;
    fwddr = A(10, :)./N;
    fwrww = A(11, :)./N;
    fwrwd = A(12, :)./N;
    fwrwr = A(13, :)./N;
    fwrdd = A(14, :)./N;
    fwrdr = A(15, :)./N;
    fddww = A(16, :)./N;
    fddwd = A(17, :)./N;
    fddwr = A(18, :)./N;
    fdddd = A(19, :)./N;
    fdddr = A(20, :)./N;
    fdrww = A(21, :)./N;
    fdrwd = A(22, :)./N;
    fdrwr = A(23, :)./N;
    fdrdd = A(24, :)./N;
    fdrdr = A(25, :)./N;
    
    %% reaction
    P = zeros(25, Nr);
    P(1, :) = fwwww.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwwwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(2, :) = 0;
    P(3, :) = fwwww.*(0.5*fwwwr+0.25*fwrwr)+fwwwr.*(0.5*fwwww+0.5*fwwwr+0.25*fwrww+0.25*fwrwr)+fwdww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrww.*(0.25*fwwwr+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.25*fwwwr+0.125*fwrww+0.125*fwrwr);
    P(4, :) = 0;
    P(5, :) = 0;
    P(6, :) = 0;
    P(7, :) = fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwddd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fddwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdddd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fdddr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrdd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(8, :) = fwdww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwdwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fddww.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fddwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fddwr.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fdddr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(9, :) = 0;
    P(10, :) = fwdwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwddd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwddr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fddwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdddd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fdddr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fdrdd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrdr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr);
    P(11, :) = fwwww.*(0.5*fwrww+0.25*fwrwr)+fwwwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwwwr.*(0.25*fwrww+0.125*fwrwr)+fwrww.*(0.5*fwwww+0.25*fwwwr+0.5*fwrww+0.25*fwrwr)+fwrwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.125*fwwwr+0.25*fwrww+0.125*fwrwr);
    P(12, :) = fwwwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwwdd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwwdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwddd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwrwd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwrdd.*(1.0*fwwww+0.5*fwwwr+0.5*fwrww+0.25*fwrwr)+fwrdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrdd.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(13, :) = fwwww.*(0.25*fwrwr)+fwwwd.*(0.25*fwwwr+0.125*fwrwr)+fwwwr.*(0.25*fwrww+0.25*fwrwr)+fwwdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fwdww.*(0.25*fwrww+0.125*fwrwr)+fwdwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwdwr.*(0.25*fwrww+0.125*fwrwr)+fwddr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fwrww.*(0.25*fwwwr+0.25*fwrwr)+fwrwd.*(0.25*fwwwr+0.125*fwrwr)+fwrwr.*(0.25*fwwww+0.25*fwwwr+0.25*fwrww+0.25*fwrwr)+fwrdr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrww.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrwd.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr)+fdrwr.*(0.5*fwwww+0.25*fwwwr+0.25*fwrww+0.125*fwrwr)+fdrdr.*(0.25*fwwww+0.125*fwwwr+0.125*fwrww+0.0625*fwrwr);
    P(14, :) = 0;
    P(15, :) = fwwwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwwdd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fwwdr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwdwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwddd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwddr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fwrwd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fwrdd.*(0.5*fwwwr+0.5*fwdww+0.5*fwdwr+0.25*fwrwr)+fwrdr.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrwd.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr)+fdrdd.*(0.25*fwwwr+0.25*fwdww+0.25*fwdwr+0.125*fwrwr)+fdrdr.*(0.125*fwwwr+0.125*fwdww+0.125*fwdwr+0.0625*fwrwr);
    P(16, :) = 0;
    P(17, :) = 0;
    P(18, :) = 0;
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
end
nd = 2*sum(A(16:20,:))+sum(A(6:10,:)+A(21:25,:));
final = nd*ring_area;
res = final / initial;
end