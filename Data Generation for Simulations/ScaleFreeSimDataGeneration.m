
%% To generate data for scale free networks of p=40 sizes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in scale free "truths" from ScaleFreeSimTruths_p40.R

Omega1_g = csvread('A1_g_p40_similar.csv');
Omega2_g = csvread('A2_g_p40_similar.csv');
Omega3_g = csvread('A3_g_p40_similar.csv');

Omega1_m = csvread('A1_m_p40_different.csv');
Omega2_m = csvread('A2_m_p40_different.csv');
Omega3_m = csvread('A3_m_p40_different.csv');

% True covariance matrix is inverse of precision matrix
Cov1_True_m = inv(Omega1_m);
Cov2_True_m = inv(Omega2_m);
Cov3_True_m = inv(Omega3_m);

% True covariance matrix is inverse of precision matrix
Cov1_True_g = inv(Omega1_g);
Cov2_True_g = inv(Omega2_g);
Cov3_True_g = inv(Omega3_g);

p_genes=40;
p_metabolites=40;
% Simulate data with sample size 100 per group
n = 100;

X1_g = rMNorm(zeros(p_genes, 1), Cov1_True_g, n)';
X2_g = rMNorm(zeros(p_genes, 1), Cov2_True_g, n)';
X3_g = rMNorm(zeros(p_genes, 1), Cov3_True_g, n)';
X1_m = rMNorm(zeros(p_metabolites, 1), Cov1_True_m, n)';
X2_m = rMNorm(zeros(p_metabolites, 1), Cov2_True_m, n)';
X3_m = rMNorm(zeros(p_metabolites, 1), Cov3_True_m, n)';

%Write out data for further simulations
csvwrite('X1_g_p40_setting1_sim1.csv', X1_g);
csvwrite('X2_g_p40_setting1_sim1.csv', X2_g);
csvwrite('X3_g_p40_setting1_sim1.csv', X3_g);
csvwrite('X1_m_p40_setting1_sim1.csv', X1_m);
csvwrite('X2_m_p40_setting1_sim1.csv', X2_m);
csvwrite('X3_m_p40_setting1_sim1.csv', X3_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in scale free truths for Setting Two
Omega1_g = csvread('A1_g_p40_different.csv');
Omega2_g = csvread('A2_g_p40_different.csv');
Omega3_g = csvread('A3_g_p40_different.csv');

Omega1_m = csvread('A1_m_p40_different.csv');
Omega2_m = csvread('A2_m_p40_different.csv');
Omega3_m = csvread('A3_m_p40_different.csv');

% True covariance matrix is inverse of precision matrix
Cov1_True_m = inv(Omega1_m);
Cov2_True_m = inv(Omega2_m);
Cov3_True_m = inv(Omega3_m);

% True covariance matrix is inverse of precision matrix
Cov1_True_g = inv(Omega1_g);
Cov2_True_g = inv(Omega2_g);
Cov3_True_g = inv(Omega3_g);

p_genes=40;
p_metabolites=40;
% Simulate data with sample size 100 per group
n = 100;

X1_g = rMNorm(zeros(p_genes, 1), Cov1_True_g, n)';
X2_g = rMNorm(zeros(p_genes, 1), Cov2_True_g, n)';
X3_g = rMNorm(zeros(p_genes, 1), Cov3_True_g, n)';
X1_m = rMNorm(zeros(p_metabolites, 1), Cov1_True_m, n)';
X2_m = rMNorm(zeros(p_metabolites, 1), Cov2_True_m, n)';
X3_m = rMNorm(zeros(p_metabolites, 1), Cov3_True_m, n)';

%Write out data for further simulations
csvwrite('X1_g_p40_setting2_sim1.csv', X1_g);
csvwrite('X2_g_p40_setting2_sim1.csv', X2_g);
csvwrite('X3_g_p40_setting2_sim1.csv', X3_g);
csvwrite('X1_m_p40_setting2_sim1.csv', X1_m);
csvwrite('X2_m_p40_setting2_sim1.csv', X2_m);
csvwrite('X3_m_p40_setting2_sim1.csv', X3_m);

% PRIOR TO ANALYSIS, DATA FILES WERE THEN STANDARDIZED TO HAVE
% A STANDARD DEVIATION OF ONE