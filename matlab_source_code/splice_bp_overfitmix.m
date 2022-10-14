function [] = splice_csp_rev_overfitmix(K, B)

rng(0513)

% load datasets_rep25_K4d4_n1000p20;

addpath('./rmvnrnd')

alpha0 = 5;

% [n,d,p,rep] = size(Yn_big_arr);

Y_with_type = csvread('slice_no_n.csv');
Y_mat_d4 = Y_with_type(:, 2:end);
gene_type = Y_with_type(:, 1);
Y_mat = Y_mat_d4;

[n, p] = size(Y_mat);
d = max(max(Y_mat));
Y_arr = zeros(n,p,d);
for c=1:d
    Y_arr(:,:,c) = (Y_mat==c);
end

% MCMC setup
nrun = 15000; burn = 5000; thin = 5; 
% nrun = 40; burn = 0; thin = 1; 

% global
% B = 3; 

% [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
%     z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
%     = sglca_csp_fun(Y_arr, K, B, nrun, alpha0);

[Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
    z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
    = sglca_csp_fun_overfitmix(Y_arr, K, B, nrun, alpha0);
% Y_arr, K, B, nrun, alpha0

K_star = [median(K_star_arr(burn+1:thin:end)); iqr(K_star_arr(burn+1:thin:end))];

z_beta_vec = mean(z_beta_vec_arr(:,burn+1:thin:end), 2);

% store other posterior means
% Q-matrix, permuted
Q_mat_pomean = mean(Q_mat_arr(:,:,burn+1:thin:end), 3);

% beta_mat, permuted
beta_mat_pomean = mean(beta_mat_arr(:,:,:,burn+1:thin:end), 4);

% get beta_tilde
beta_mat_tilde_arr = zeros(p, K, d, nrun);
for ii=1:nrun
    beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
end
beta_mattil_pomean = mean(beta_mat_tilde_arr(:,:,:,burn+1:thin:end), 4);

% beta0 posterior mean
beta0_pomean = mean(beta0_arr(:,:,burn+1:thin:end), 3);

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);

% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:thin:end), 2);

% A_mat and z posterior mean; A_mat permuted
A_mat_pomean = mean(A_mat_arr(:,:,burn+1:thin:end), 3);
z_mat_pomean = mean(z_tau_mat_arr(:,:,burn+1:thin:end), 3);
% permute_vec_rep(g,:) = permute_vec;
sig2_beta = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);


filename = strcat('s1plice_overfitmix_n', num2str(n), 'p', num2str(p),...
    '_Kini', num2str(K), '_B', num2str(B), '_d', num2str(d), '.mat');

save(filename);

end
