function [] = bp_csp_simu(K, alpha0)

% DEBUG mode
% K = 7; alpha0 = 5; g = 5;

rng(0513)

% set up parpool
% thePool = parpool(25);

% add path to a truncated normal sampler
addpath('./rmvnrnd')

filename = 'datasets_strong_rep50_K4d4_n2000p20.mat';
load(filename);

[n,d,p,rep] = size(Yn_big_arr);

% MCMC setup
nrun = 15000; burn = 5000; thin = 5; 

% global
B = 2;

% store for replications
Q_mat_pomean_rep = zeros(p,K,rep);
%
beta0_pomean_rep = zeros(p,d,rep);
beta_mat_pomean_rep = zeros(p,K,d,rep);
beta_mattil_pomean_rep = zeros(p,K,d,rep);
%
Bern_K_pomean_rep = zeros(K,B,rep);
tau_pomean_rep = zeros(B,rep);
%
A_mat_pomean_rep = zeros(n,K,rep);
z_mat_pomean_rep = zeros(n,B,rep);
%
DIC_rep = zeros(1,nrun);
% permute_vec_rep = zeros(rep,K);
sig2_beta_rep = zeros(K,d-1,rep);
% CSP-specific quantities
K_star_rep = zeros(2,rep);
z_beta_vec_rep = zeros(K,rep);

% NEW: track K_star_arr
K_star_arrs = zeros(rep, nrun);

% NEW: debug
beta_mat_arrs = cell(rep,1);
beta0_arrs = cell(rep,1);
Q_mat_arrs = cell(rep,1);
sig2_beta_arrs = cell(rep,1);

% Use parallel computing to execute independent simulation replicates
for g = 1:rep
            
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

    % -- initialize parameters -- %
    [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
        z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
        = bp_csp_fun(Y_arr, K, B, nrun, alpha0);
    
    K_star_rep(:, g) = [median(K_star_arr(burn+1:thin:end)); iqr(K_star_arr(burn+1:thin:end))];
    % track K_star_arr
    K_star_arrs(g, :) = K_star_arr;
    
    z_beta_vec_rep(:, g) = mean(z_beta_vec_arr(:,burn+1:thin:end), 2);
    
    % store other posterior means
    % Q-matrix, permuted
    Q_mat_pomean_rep(:,:,g) = mean(Q_mat_arr(:,:,burn+1:thin:end), 3);
    
    % beta_mat, permuted
    beta_mat_pomean_rep(:,:,:,g) = mean(beta_mat_arr(:,:,:,burn+1:thin:end), 4);
    
    % get beta_tilde
    beta_mat_tilde_arr = zeros(p, K, d, nrun);
    for ii=1:nrun
        beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
    end
    beta_mattil_pomean_rep(:,:,:,g) = mean(beta_mat_tilde_arr(:,:,:,burn+1:thin:end), 4);
    
    % beta0 posterior mean
    beta0_pomean_rep(:,:,g) = mean(beta0_arr(:,:,burn+1:thin:end), 3);
    
    % Bern_K posterior mean
    Bern_K_pomean_rep(:,:,g) = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);

    % tau posterior mean
    tau_pomean_rep(:,g) = mean(tau_arr(:,burn+1:thin:end), 2);
    
    % A_mat and z posterior mean; A_mat permuted
    A_mat_pomean_rep(:,:,g) = mean(A_mat_arr(:,:,burn+1:thin:end), 3);
    z_mat_pomean_rep(:,:,g) = mean(z_tau_mat_arr(:,:,burn+1:thin:end), 3);
    sig2_beta_rep(:,:,g) = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
        
    % NEW
    beta_mat_arrs{g} = beta_mat_arr;
    beta0_arrs{g} = beta0_arr;
    Q_mat_arrs{g} = Q_mat_arr;
    sig2_beta_arrs{g} = sig2_beta_arr;

end


filename = strcat('msimu_CSPbeta_n', num2str(n), 'p', num2str(p),...
    '_Kini', num2str(K), '_d', num2str(d), '_rep', num2str(rep), '.mat');

save(filename);


% delete(thePool);


end