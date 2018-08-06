addpath './Wang Li code - original download';

%Example: Call MCMC sampler for multiple graphs and multiple platforms
% on simulated data with K=3 subgroups and S=2 platforms where platform 1
% has 15 variables and platform 2 has 20 variables

% p is number of variables
p_1 = 15;
p_2 = 20;
p = [p_1 p_2];

% K is number of subgroups
K = 3;

% S is number of platforms
S = 2;


% Set up a precision matrix for each group for platform 1.
% Here they are all the same.
Omega1_1 = toeplitz([1, 0.5, 0.4, zeros(1, p_1 - 3)]);
Omega2_1 = Omega1_1;
Omega3_1 = Omega1_1;

% True covariance matrix is inverse of precision matrix
Cov1_True_1 = inv(Omega1_1);
Cov2_True_1 = inv(Omega2_1);
Cov3_True_1 = inv(Omega3_1);

%Set up precision matrices for platform 2
%Here they are all different
Omega1_2 = toeplitz([1,0.5,zeros(1,p_2-2)]);
perm_group2 = randperm(p_2);
perm_group3 = randperm(p_2);
Omega2_2 = Omega1_2(perm_group2, perm_group2);
Omega3_2 = Omega1_2(perm_group3, perm_group3);

% True covariance matrix is inverse of precision matrix
Cov1_True_2 = inv(Omega1_2);
Cov2_True_2 = inv(Omega2_2);
Cov3_True_2 = inv(Omega3_2);

% Save out these matrices to use as input to simulation script
%csvwrite('A1_1.csv', Omega1_1)
%csvwrite('A2_1.csv', Omega2_1)
%csvwrite('A3_1.csv', Omega3_1)
%csvwrite('A1_2.csv', Omega1_2)
%csvwrite('A2_2.csv', Omega2_2)
%csvwrite('A3_2.csv', Omega3_2)

% Simulate data with sample size 100 per group
n = 100;

% Random normal data with mean 0 and given covariance and sample size
X1_1 = rMNorm(zeros(p_1, 1), Cov1_True_1, n)';
X2_1 = rMNorm(zeros(p_1, 1), Cov2_True_1, n)';
X3_1 = rMNorm(zeros(p_1, 1), Cov3_True_1, n)';
X1_2 = rMNorm(zeros(p_2, 1), Cov1_True_2, n)';
X2_2 = rMNorm(zeros(p_2, 1), Cov2_True_2, n)';
X3_2 = rMNorm(zeros(p_2, 1), Cov3_True_2, n)';

%Write out data for further simulations
%csvwrite('X1_1.csv', X1_1)
%csvwrite('X2_1.csv', X2_1)
%csvwrite('X3_1.csv', X3_1)
%csvwrite('X1_2.csv', X1_2)
%csvwrite('X2_2.csv', X2_2)
%csvwrite('X3_2.csv', X3_2)

%X1_1 = csvread('X1_1.csv');
%X2_1 = csvread('X2_1.csv');
%X3_1 = csvread('X3_1.csv');
%X1_2 = csvread('X1_2.csv');
%X2_2 = csvread('X2_2.csv');
%X3_2 = csvread('X3_2.csv');

% X'X matrix
S1_1 = X1_1' * X1_1;
S2_1 = X2_1' * X2_1;
S3_1 = X3_1' * X3_1;
S1_2 = X1_2' * X1_2;
S2_2 = X2_2' * X2_2;
S3_2 = X3_2' * X3_2;

SampCov_1 = cat(3, S1_1, S2_1, S3_1);
SampCov_2 = cat(3, S1_2, S2_2, S3_2);
SampCov = {SampCov_1; SampCov_2};

% Number of MCMC iterations before burnin
burnin  = 10000;

% Number of MCMC iterations after burnin
nmc = 10000;

%% Setting Hyperparameters %%

% Increasing h below results in more sparsity of the inferred network
h = 50^2;
v0=0.02^2;
v1=h*v0;
lambda = 1;
pii_1 = 2/(p_1-1);
pii_2 = 2/(p_2-1);
pii_1 = repmat(pii_1, 1, K);
pii_2 = repmat(pii_2, 1, K);
pii = {pii_1; pii_2};

% Graph Specific Hyperparameters %
V0_vector_1 = v0*ones(p_1);
V1_vector_1 = v1*ones(p_1);
V0_1 = cat(3, V0_vector_1, V0_vector_1, V0_vector_1);
V1_1 = cat(3, V1_vector_1, V1_vector_1, V1_vector_1);

V0_vector_2 = v0*ones(p_2);
V1_vector_2 = v1*ones(p_2);
V0_2 = cat(3, V0_vector_2, V0_vector_2, V0_vector_2);
V1_2 = cat(3, V1_vector_2, V1_vector_2, V1_vector_2);

V0 = {V0_1; V0_2};
V1 = {V1_1; V1_2};

% Network accuracy is stable across varying hyperparameter settings,
% however relative similarity measures Theta and Phi may vary in 
% magnitude (order remains fairly consistent) as alpha, beta, eta, kappa,
% a, b, d, f change
%
% For gamma slab priors, alpha & beta, eta & kappa, appropriate choices are
% those which allow good discrimination between zero and non-zero values of
% Theta (determined by alpha and beta) and Phi (determined by eta and
% kappa). A possible prior is one of the form Gamma(x|alpha, beta) with 
% alpha>1 as this prior is equal to 0 at x=0 and non-zero for x>0 
%
% Prior parameters for gamma slab of mixture prior for groups
alpha = repmat(1, 1, S);
beta = 9;

%Prior parameters for gamma slab of mixture prior for platforms
eta = 4;
kappa = 5;

% Recommendations for setting a and b are to set a and b such that the
% resulting prior probability of edge inclusion is smaller than the
% expected level of sparsity

% Parameters for prior on nu which affects graph sparsity for groups
a = 1;
b = 7;

% Parameters for prior on w which affects graph sparsity for platform
% relatedness
d = 1;
f = 19;

% Initial value for precision matrix for each group
InitialC_1 = eye(p_1);
InitialC_2 = eye(p_2);

% Initial values for Theta and nu
Theta = zeros(K);
Theta = {Theta; Theta};

nu_1 = zeros(p_1, p_1) - 1;
nu_2 = zeros(p_2, p_2) - 1;
nu = {nu_1; nu_2};

% Initial values for precision and covariance matrices 
C_1 = cat(3, InitialC_1, InitialC_1, InitialC_1);
C_2 = cat(3, InitialC_2, InitialC_2, InitialC_2);
C = {C_1; C_2};

% When accessing cells from cell array
% C{1,:} returns C_1
% C{1,:}(:,:,1) returns the first groups initial matrix for platform 1

InitialSig_1 = inv(InitialC_1);
Sig_1 = cat(3, InitialSig_1, InitialSig_1, InitialSig_1);
InitialSig_2 = inv(InitialC_2);
Sig_2 = cat(3, InitialSig_2, InitialSig_2, InitialSig_2);
Sig = {Sig_1; Sig_2};

% Initial values for Phi and w
Phi = zeros(S);
my_w = zeros(K,K)-1;

% Number of Subjects per platform
n_1 = repmat(n, 1, K);
n_2 = repmat(n, 1, K);
n = {n_1; n_2};

% Parameter for Bernoulli prior on indicator of platform graph relatedness
my_u = .1;


% Multiple Platform Test
[C_save, Sig_save, adj_save, Phi_save, Theta_save, ar_gamma, ar_theta, ar_zeta, ar_phi, nu_save, ar_nu, w_save, ar_w] = ...
    MCMC_multiple_graphs_SSVS_with_platforms(Theta, Phi, Sig, V0, V1, lambda, pii, n, SampCov, ...
    C, nu, alpha, beta, a, b, my_w, d, f, eta, kappa, my_u, burnin, nmc, true);


%PPIs for Phi (platform similarity measure)
ppi_phi = mean(Phi_save ~= 0, 3);

% PPIs for Theta (graph similarity measure)
ppi_theta_1 = mean(Theta_save{1,:} ~= 0, 3);
ppi_theta_2 = mean(Theta_save{2,:} ~= 0, 3);

% Edge PPIs for each graph
ppi_edges_1 = mean(adj_save{1,:}, 4);
ppi_edges_2 = mean(adj_save{2,:}, 4);

% Adjacency Matrices giving indication of edge locations for each graph
edges_1 = mean(adj_save{1,:}, 4)>.5;
edges_2 = mean(adj_save{2,:}, 4)>.5;

% Get 95% credible intervals for omega (precision matrix)
CI_omega_lower_1 = quantile(C_save{1}, 0.025, 4);
CI_omega_upper_1 = quantile(C_save{1}, 0.975, 4);
CI_omega_lower_2 = quantile(C_save{2}, 0.025, 4);
CI_omega_upper_2 = quantile(C_save{2}, 0.975, 4);