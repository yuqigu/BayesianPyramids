function [] = sglca_simu(K)

rng(0513)

load datasets_rep25_K4d4_n1000p20;

addpath('./rmvnrnd')

[n,d,p,rep] = size(Yn_big_arr); 

% global
B = 2; 

% parpool setup
parpool(rep);
pctRunOnAll load('datasets_rep25_K4d4_n1000p20.mat')


% MCMC setup
nrun = 15000; burn = 5000; thin = 5; 

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
% permute_vec_rep = zeros(rep,K);
sig2_beta_rep = zeros(K,d-1,rep);


DIC_rep = zeros(rep,1);

parfor g=1:rep
        
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

    % -- initialize parameters -- %
    [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr, z_mat_arr, sig2_beta_arr] ...
        = sglca_fun(Y_arr, K, B, nrun, burn, thin);
    
    % store other posterior means
    % Q-matrix, permuted
    % Q_mat_pomean_rep(:,:,g) = Q_pm_raw(:, permute_vec);
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
    z_mat_pomean_rep(:,:,g) = mean(z_mat_arr(:,:,burn+1:thin:end), 3);
    % permute_vec_rep(g,:) = permute_vec;
    sig2_beta_rep(:,:,g) = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
    
    %% compute the DIC value
    eff_nrun = (nrun-burn)/thin;
    loglik_data_eff = zeros(n, eff_nrun);
    penalty_eff = zeros(n, eff_nrun);
    for ii=burn+1:thin:nrun
        ind = ceil((ii-burn)/thin);
        loglik_data_eff(:,ind) = loglik(Y_arr, Q_mat_arr(:,:,ind), beta_mat_arr(:,:,:,ind), beta0_arr(:,:,ind), ...
            Bern_K_arr(:,:,ind), tau_arr(:,ind), A_mat_arr(:,:,ind), z_mat_arr(:,:,ind));
        
        Q_rep_g = (Q_mat_pomean_rep(:,:,g) > 0.5);
        penalty_eff(:,ind) = loglik(Y_arr, Q_rep_g, beta_mat_pomean_rep(:,:,:,g), beta0_pomean_rep(:,:,g), ...
            Bern_K_pomean_rep(:,:,g), tau_pomean_rep(:,g), A_mat_arr(:,:,ind), z_mat_arr(:,:,ind));
    end
    
    DIC_rep(g) = - 4/eff_nrun * sum(sum(loglik_data_eff)) ...
        + 2/eff_nrun * sum(sum(penalty_eff));

end

filename = strcat('msimu_DICsig_n', num2str(n), 'p', num2str(p),...
    '_Kini', num2str(K), '_d', num2str(d), '_rep', num2str(rep), '.mat');

save(filename);

end




