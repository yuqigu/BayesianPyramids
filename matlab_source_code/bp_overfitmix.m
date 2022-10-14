n = 1000;
filename = strcat('datasets_rep50_K4d4_n', num2str(n), 'p20.mat');
load(filename);


% add path to a truncated normal sampler
addpath('./rmvnrnd')

[n,d,p,rep] = size(Yn_big_arr);


% MCMC setup
nrun = 15000; burn = 5000; thin = 5;
% nrun = 1500; burn = 1000; thin = 1; 


% global
B = 5; K = 4;
alpha0 = 5;

% store for replications
Q_mat_pomean_rep = zeros(p,K,rep);
%
A_mat_pomean_rep = zeros(n,K,rep);
z_mat_pomean_rep = zeros(n,B,rep);
%
% DIC_rep = zeros(1,nrun);
% permute_vec_rep = zeros(rep,K);
sig2_beta_rep = zeros(K,d-1,rep);
% CSP-specific quantities
K_star_rep = zeros(2,rep);
z_beta_vec_rep = zeros(K,rep);

% track K_star_arr
K_star_arrs = zeros(rep, nrun);

%
beta0_pomean_rep = zeros(p,d,rep);
beta_mat_pomean_rep = zeros(p,K,d,rep);
beta_mattil_pomean_rep = zeros(p,K,d,rep);
%
Bern_K_pomean_rep = zeros(K,B,rep);
tau_pomean_rep = zeros(B,rep);

L = (nrun - burn)/thin;

tic;

% set up parpool
% parpool(20);
parfor g = 1:rep
% for g = 1:1
            
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

    [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
        z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
        = sglca_csp_fun_overfitmix(Y_arr, K, B, nrun, alpha0);
    
    K_star_rep(:, g) = [median(K_star_arr(burn+1:thin:end)); iqr(K_star_arr(burn+1:thin:end))];
    
    z_beta_vec_rep(:, g) = mean(z_beta_vec_arr(:,burn+1:thin:end), 2);
    
    % store other posterior means
    % Q-matrix, permuted
    Q_mat_pomean_rep(:,:,g) = mean(Q_mat_arr(:,:,burn+1:thin:end), 3);
    
    % A_mat and z posterior mean; A_mat permuted
    A_mat_pomean_rep(:,:,g) = mean(A_mat_arr(:,:,burn+1:thin:end), 3);
    z_mat_pomean_rep(:,:,g) = mean(z_tau_mat_arr(:,:,burn+1:thin:end), 3);
    % permute_vec_rep(g,:) = permute_vec;
    sig2_beta_rep(:,:,g) = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
    
    
    % -- continuous parameters begin here -- %
    % beta_mat, permuted
    beta_mat_pomean_rep(:,:,:,g) = mean(beta_mat_arr(:,:,:,burn+1:thin:end), 4);
    
    % get beta_tilde
    beta_mat_tilde_arr = zeros(p, K, d, nrun);
    for ii=1:nrun
        beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
    end
    
    % beta_mat posterior mean
    beta_mattil_pomean_rep(:,:,:,g) = mean(beta_mat_tilde_arr(:,:,:,burn+1:thin:end), 4);
    
    % beta0 posterior mean
    beta0_pomean_rep(:,:,g) = mean(beta0_arr(:,:,burn+1:thin:end), 3);
    
    % Bern_K posterior mean
    Bern_K_pomean_rep(:,:,g) = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);

    % tau posterior mean
    tau_pomean_rep(:,g) = mean(tau_arr(:,burn+1:thin:end), 2);
    
    
    fprintf('Replicate %d completed\n', g);
end

save('mdata_rev_bp_overfitmix.mat')

mytime = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_pomean_sort = sort(tau_pomean_rep, 1, 'descend');

figure;
boxplot(tau_pomean_sort')
set(gca, 'FontSize', 16)
pbaspect([3 2 1])
title('Deep latent class proportions across simulation replicates')
xlabel('B=5 deep latent classes')
ylabel('posterior means of proportions')
% print('-r300', 'overfitmix_simu', '-dpng')



