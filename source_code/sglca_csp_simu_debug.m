function [] = sglca_csp_simu_debug(K, alpha0)

% How about we use CSP prior for the variance of beta's? 
% only put CSP prior on beta_{j1k} for c=1

% DEBUG mode
% K = 7; alpha0 = 5; g = 5;

rng(0513)

% % -- for greatlakes begin -- %
% %%%%  We get from the environment the number of processors
% NP = str2num(getenv('SLURM_NTASKS'));
% % parpool setup
% %%%%  Create the pool for parfor to use
% thePool = parpool('current', NP);
% % -- for greatlakes end -- %

load datasets_rep50_K4d4_n1000p20;

% % set up parpool
% thePool = parpool(36);
% pctRunOnAll load('datasets_rep50_K4d4_n1000p20.mat')


% add path to a truncated normal sampler
addpath('./rmvnrnd')


[n,d,p,rep] = size(Yn_big_arr);


% MCMC setup
% nrun = 15000; burn = 5000; thin = 5; 
nrun = 500; burn = 100; thin = 1; 


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


% parfor g=1:rep
        
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

    % -- initialize parameters -- %
%     [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr, z_tau_mat_arr, sig2_beta_arr] ...
%         = sglca_fun(Y_arr, K, B, nrun, burn, thin);

    [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
        z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
        = sglca_csp_fun_allpos(Y_arr, K, B, nrun, alpha0);
    
    K_star_rep(:, g) = [median(K_star_arr(burn+1:thin:end)); iqr(K_star_arr(burn+1:thin:end))];
    
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
    % permute_vec_rep(g,:) = permute_vec;
    sig2_beta_rep(:,:,g) = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
    
    
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
    
    DIC_rep(g) = - 4/eff_nrun * sum(sum(loglik_data_eff)) ...
        + 2/eff_nrun * sum(sum(penalty_eff));

% end


% filename = strcat('local_msimu_CSPbeta_n', num2str(n), 'p', num2str(p),...
%     '_Kini', num2str(K), '_d', num2str(d), '_rep', num2str(rep), '.mat');
% 
% save(filename);


% delete(thePool);


end