addpath '/Users/elinsellers/Desktop/Multiple graphs SSVS code/Wang Li code - original download';


%% To set up simulation settings for AR(2) networks of p=80 sizes
p_genes = 80;
p_metabolites = 80;
% K is number of subgroups
K = 3;

% S is number of platforms
S = 2;

% Set up a precision matrix for each group for genes. Here they are all
% similar, with about 10% of edges varying between group 1 and 2, and group
% 1 and 3
Omega1_g = toeplitz([1, 0.5, 0.4, zeros(1, p_genes - 3)]);
n_edges = (sum(sum(Omega1_g ~= 0)) - p_genes) / 2;  %157
n_possible = p_genes * (p_genes - 1) / 2; %3160
change = randsample(p_genes,6) %30, 40 | 54, 24, 79, 6 | 69, 33, 27, 51, 8 | 3 61 69 19 57 5
perm_2=1:80; %sum is 3240
perm_2(8)=69;
perm_2(27)=33;
perm_2(33)=27;
perm_2(51)=51;
perm_2(69)=8;
perm_3=1:80;
perm_3(3)=3;
perm_3(5)=61;
perm_3(19)=69;
perm_3(57)=19;
perm_3(61)=57;
perm_3(69)=5;
Omega2_g = Omega1_g(perm_2, perm_2);
Omega3_g = Omega1_g(perm_3, perm_3);
adj1_g = Omega1_g>.1;
adj2_g = Omega2_g>.1;
adj3_g = Omega3_g>.1;
count_shared_edges(adj1_g,adj2_g) %141
count_shared_edges(adj1_g, adj3_g) %138
count_shared_edges(adj2_g, adj3_g) %128


% Set up a precision matrix for each group for metabolites. Here they are all
% similar, with about 10% of edges varying between group 1 and 2, and group
% 1 and 3
Omega1_m = toeplitz([1, 0.5, 0.4, zeros(1, p_metabolites - 3)]);
n_edges = (sum(sum(Omega1_m ~= 0)) - p_genes) / 2;  %157
n_possible = p_metabolites * (p_metabolites - 1) / 2; %3160
change = randsample(p_metabolites,5) % 38, 42, 18, 13, 41
change = randsample(p_metabolites,5) % 60, 73, 21, 55, 45
perm_2=1:80; %sum is 3240
perm_2(13)=38;
perm_2(18)=42;
perm_2(38)=18;
perm_2(41)=13;
perm_2(42)=41;
perm_3=1:80;
perm_3(21)=60;
perm_3(45)=73;
perm_3(55)=21;
perm_3(60)=55;
perm_3(73)=45;
Omega2_m = Omega1_m(perm_2, perm_2);
Omega3_m = Omega1_m(perm_3, perm_3);
adj1_m = Omega1_m>.1;
adj2_m = Omega2_m>.1;
adj3_m = Omega3_m>.1;
count_shared_edges(adj1_m,adj2_m) %140
count_shared_edges(adj1_m, adj3_m) %137
count_shared_edges(adj2_m, adj3_m) %120


% Save out these matrices to use as input to simulation script for BOTH
% PLATFORMS SIMILAR
csvwrite('A1_g_p80_similar.csv', Omega1_g)
csvwrite('A2_g_p80_similar.csv', Omega2_g)
csvwrite('A3_g_p80_similar.csv', Omega3_g)
csvwrite('A1_m_p80_similar.csv', Omega1_m)
csvwrite('A2_m_p80_similar.csv', Omega2_m)
csvwrite('A3_m_p80_similar.csv', Omega3_m)

imagesc(adj1_m)
imagesc(adj2_m)
imagesc(adj3_m)
imagesc(adj1_g)
imagesc(adj2_g)
imagesc(adj3_g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up a precision matrix for each group for genes. Here they are all
% different
Omega1_g = toeplitz([1, 0.5, 0.4, zeros(1, p_genes - 3)]);
perm_group2 = randperm(p_genes);
perm_group3 = randperm(p_genes);
Omega2_g = Omega1_g(perm_group2, perm_group2);
Omega3_g = Omega1_g(perm_group3, perm_group3);
adj1_g = Omega1_g>.1;
adj2_g = Omega2_g>.1;
adj3_g = Omega3_g>.1;
count_shared_edges(adj1_g,adj2_g) %8
count_shared_edges(adj1_g, adj3_g) %6
count_shared_edges(adj2_g, adj3_g) %6

% Set up a precision matrix for each group for metabolites. Here they are all
% different
Omega1_m = toeplitz([1, 0.5, 0.4, zeros(1, p_metabolites - 3)]);
perm_group2 = randperm(p_metabolites);
perm_group3 = randperm(p_metabolites);
Omega2_m = Omega1_m(perm_group2, perm_group2);
Omega3_m = Omega1_m(perm_group3, perm_group3);
adj1_m = Omega1_m>.1;
adj2_m = Omega2_m>.1;
adj3_m = Omega3_m>.1;
count_shared_edges(adj1_m,adj2_m) %8
count_shared_edges(adj1_m, adj3_m) %6
count_shared_edges(adj2_m, adj3_m) %4


% Save out these matrices to use as input to simulation script for BOTH
% PLATFORMS DIFFERENT
csvwrite('A1_g_p80_different.csv', Omega1_g)
csvwrite('A2_g_p80_different.csv', Omega2_g)
csvwrite('A3_g_p80_different.csv', Omega3_g)
csvwrite('A1_m_p80_different.csv', Omega1_m)
csvwrite('A2_m_p80_different.csv', Omega2_m)
csvwrite('A3_m_p80_different.csv', Omega3_m)

imagesc(adj1_m)
imagesc(adj2_m)
imagesc(adj3_m)
imagesc(adj1_g)
imagesc(adj2_g)
imagesc(adj3_g)

% To check matrix is positive definite
% If the input matrix is not positive definite, then the output will be a positive integer
[~,pos]=chol(Omega1_m) %0
[~,pos]=chol(Omega2_m) %0
[~,pos]=chol(Omega3_m) %0
[~,pos]=chol(Omega1_g) %0
[~,pos]=chol(Omega2_g) %0
[~,pos]=chol(Omega3_g) %0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To generate sample test case for setting one
Omega1_g = csvread('A1_g_p80_similar.csv');
Omega2_g = csvread('A2_g_p80_similar.csv');
Omega3_g = csvread('A3_g_p80_similar.csv');
Omega1_m = csvread('A1_m_p80_different.csv');
Omega2_m = csvread('A2_m_p80_different.csv');
Omega3_m = csvread('A3_m_p80_different.csv');

% True covariance matrix is inverse of precision matrix
Cov1_True_m = inv(Omega1_m);
Cov2_True_m = inv(Omega2_m);
Cov3_True_m = inv(Omega3_m);

% True covariance matrix is inverse of precision matrix
Cov1_True_g = inv(Omega1_g);
Cov2_True_g = inv(Omega2_g);
Cov3_True_g = inv(Omega3_g);

% Simulate data with sample size 100 per group
n = 100;

X1_g = rMNorm(zeros(p_genes, 1), Cov1_True_g, n)';
X2_g = rMNorm(zeros(p_genes, 1), Cov2_True_g, n)';
X3_g = rMNorm(zeros(p_genes, 1), Cov3_True_g, n)';
X1_m = rMNorm(zeros(p_metabolites, 1), Cov1_True_m, n)';
X2_m = rMNorm(zeros(p_metabolites, 1), Cov2_True_m, n)';
X3_m = rMNorm(zeros(p_metabolites, 1), Cov3_True_m, n)';


%Write out data for further simulations
csvwrite('X1_g_p80_setting1_sim1.csv', X1_g);
csvwrite('X2_g_p80_setting1_sim1.csv', X2_g);
csvwrite('X3_g_p80_setting1_sim1.csv', X3_g);
csvwrite('X1_m_p80_setting1_sim1.csv', X1_m);
csvwrite('X2_m_p80_setting1_sim1.csv', X2_m);
csvwrite('X3_m_p80_setting1_sim1.csv', X3_m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To generate test case for setting two

Omega1_m = csvread('A1_m_p80_different.csv');
Omega2_m = csvread('A2_m_p80_different.csv');
Omega3_m = csvread('A3_m_p80_different.csv');
Omega1_g = csvread('A1_g_p80_different.csv');
Omega2_g = csvread('A2_g_p80_different.csv');
Omega3_g = csvread('A3_g_p80_different.csv');

% True covariance matrix is inverse of precision matrix
Cov1_True_m = inv(Omega1_m);
Cov2_True_m = inv(Omega2_m);
Cov3_True_m = inv(Omega3_m);

% True covariance matrix is inverse of precision matrix
Cov1_True_g = inv(Omega1_g);
Cov2_True_g = inv(Omega2_g);
Cov3_True_g = inv(Omega3_g);

p_genes=80;
p_metabolites=80;
% Simulate data with sample size 100 per group
n = 100;

X1_g = rMNorm(zeros(p_genes, 1), Cov1_True_g, n)';
X2_g = rMNorm(zeros(p_genes, 1), Cov2_True_g, n)';
X3_g = rMNorm(zeros(p_genes, 1), Cov3_True_g, n)';
X1_m = rMNorm(zeros(p_metabolites, 1), Cov1_True_m, n)';
X2_m = rMNorm(zeros(p_metabolites, 1), Cov2_True_m, n)';
X3_m = rMNorm(zeros(p_metabolites, 1), Cov3_True_m, n)';


%Write out data for further simulations
csvwrite('X1_g_p80_setting2_sim1.csv', X1_g);
csvwrite('X2_g_p80_setting2_sim1.csv', X2_g);
csvwrite('X3_g_p80_setting2_sim1.csv', X3_g);
csvwrite('X1_m_p80_setting2_sim1.csv', X1_m);
csvwrite('X2_m_p80_setting2_sim1.csv', X2_m);
csvwrite('X3_m_p80_setting2_sim1.csv', X3_m);

% PRIOR TO ANALYSIS, DATA FILES WERE THEN STANDARDIZED TO HAVE
% A STANDARD DEVIATION OF ONE
