function [] = csp_splice(K)

rng(0513)

% load datasets_rep25_K4d4_n1000p20;

addpath('./rmvnrnd')

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
% nrun = 15000; burn = 5000; thin = 5; 
nrun = 200; burn = 40; thin = 1; 

% global
B = 2; 


% % store for replications
% Q_mat_pomean_rep = zeros(p,K,g);
% %
% beta0_pomean_rep = zeros(p,d,g);
% beta_mat_pomean_rep = zeros(p,K,d,g);
% beta_mattil_pomean_rep = zeros(p,K,d,g);
% %
% Bern_K_pomean_rep = zeros(K,B,g);
% tau_pomean_rep = zeros(B,g);
% %
% A_mat_pomean_rep = zeros(n,K,g);
% z_mat_pomean_rep = zeros(n,B,g);
% %
% DIC_rep = zeros(1,nrun);
% % permute_vec_rep = zeros(g,K);
% sig2_beta_rep = zeros(K,d-1,g);
% % CSP-specific quantities
% K_star_rep = zeros(2,g);
% z_beta_vec_rep = zeros(K,g);

        
% Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

% -- initialize parameters -- %
%     [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr, z_tau_mat_arr, sig2_beta_arr] ...
%         = sglca_fun(Y_arr, K, B, nrun, burn, thin);

[Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
    z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
    = sglca_csp_fun(Y_arr, K, B, nrun);

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


%% compute the DIC value
eff_nrun = (nrun-burn)/thin;
loglik_data_eff = zeros(n, eff_nrun);
penalty_eff = zeros(n, eff_nrun);
for ii=burn+1:thin:nrun
    ind = ceil((ii-burn)/thin);
    loglik_data_eff(:,ind) = loglik(Y_arr, Q_mat_arr(:,:,ind), beta_mat_arr(:,:,:,ind), beta0_arr(:,:,ind), ...
        Bern_K_arr(:,:,ind), tau_arr(:,ind), A_mat_arr(:,:,ind), z_tau_mat_arr(:,:,ind));

    Q_rep_g = (Q_mat_pomean_rep(:,:,g) > 0.5);
    penalty_eff(:,ind) = loglik(Y_arr, Q_rep_g, beta_mat_pomean_rep(:,:,:,g), beta0_pomean_rep(:,:,g), ...
        Bern_K_pomean_rep(:,:,g), tau_pomean_rep(:,g), A_mat_arr(:,:,ind), z_tau_mat_arr(:,:,ind));
end

DIC_rep = - 4/eff_nrun * sum(sum(loglik_data_eff)) ...
    + 2/eff_nrun * sum(sum(penalty_eff));


filename = strcat('splice_CSPbeta_n', num2str(n), 'p', num2str(p),...
    '_Kini', num2str(K), '_d', num2str(d), '.mat');

save(filename);

end
