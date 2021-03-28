function [] = cv_simu(K)

% performs cross validation to select the model dimension K
rng(0513)

% load datasets_rep25_K4d4_n2000p20;
load datasets_rep25_K4d4_n1000p20;

addpath('./rmvnrnd')

[n,d,p,rep] = size(Yn_big_arr); 

% % parpool setup
% parpool(rep);
% pctRunOnAll load('datasets_rep25_K4d4_n1000p20.mat')


% MCMC setup
nrun = 12000; burn = 4000; thin = 5; 

% global
B = 2; 

% % prior for intercept
% mu_beta0 = -3 * ones(3,1);
% var_beta0 = 10 * ones(d-1,1);
% 
% % prior for beta_mat
% mu_beta = zeros(K, d-1);
% 
% 
% % inverse Gamma prior for hypervariances
% a_sig = 2; b_sig = 2;
% 
% % sig2_pseudo = v0^2;
% sig2_pseudo = 0.07;


% % store for replications
% Q_mat_rep = zeros(p,K,rep);
% %
% beta0_pomean_rep = zeros(p,d,rep);
% beta_mat_pomean_rep = zeros(p,K,d,rep);
% beta_mattil_pomean_rep = zeros(p,K,d,rep);
% %
% Bern_K_pomean_rep = zeros(K,B,rep);
% tau_pomean_rep = zeros(B,rep);
% %
% A_mat_pomean_rep = zeros(n,K,rep);
% z_mat_pomean_rep = zeros(n,B,rep);
% %
% permute_vec_rep = zeros(rep,K);
% %
% DIC_rep = zeros(rep,1);


% cross validation setup
nfold = 5;
% prop_train = 1 - 1/nfold;

ll_pred_cv_rep = zeros(rep, 1);

% simulation replications begin here
for g=1:rep
    
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);
    
    % generate cross validation indices
    indices = crossvalind('Kfold', 1:n, nfold);
    
    ll_pred_vec = zeros(nfold, 1);
    
    for f=1:nfold
        Y_train_arr = Y_arr(indices ~= f, :, :); 
        Y_test_arr = Y_arr(indices == f, :, :);
    
        
        [Q_mat_arr, beta_mat_arr, beta0_arr, ~, Bern_K_arr, tau_arr, ~, ~] ...
            = sglca_fun(Y_train_arr, K, B, nrun, burn, thin);

        ll_pred_vec(f) = ll_predict(Y_test_arr, Q_mat_arr, beta_mat_arr, beta0_arr, ...
            Bern_K_arr, tau_arr, nrun, burn, thin);
    end
    ll_pred_cv_rep(g) = sum(ll_pred_vec);
    
end



filename = strcat('msimu_CV_n', num2str(n), 'p', num2str(p),...
    '_Kini', num2str(K), '_d', num2str(d), '_rep', num2str(rep), '.mat');

save(filename);


end
