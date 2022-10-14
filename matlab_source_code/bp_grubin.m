n = 1000;
filename = strcat('datasets_rep50_K4d4_n', num2str(n), 'p20.mat');
load(filename);


% add path to a truncated normal sampler
addpath('./rmvnrnd')

[n,d,p,rep] = size(Yn_big_arr);


% MCMC setup
% nrun = 15000; burn = 5000; thin = 5;
nrun = 500; burn = 300; thin = 1;


% global
B = 2; K = 7;
K0 = 4;
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

% % debug
% beta_mat_arrs = cell(rep,1);
% beta0_arrs = cell(rep,1);
% Q_mat_arrs = cell(rep,1);
% sig2_beta_arrs = cell(rep,1);


% Gelman-Rubin diagnostic for continuous parameters
nchain = 5;
%
beta0_pomean_rep = zeros(p,d,nchain,rep);
beta_mat_pomean_rep = zeros(p,K,d,nchain,rep);
beta_mattil_pomean_rep = zeros(p,K,d,nchain,rep);
%
Bern_K_pomean_rep = zeros(K,B,nchain,rep);
tau_pomean_rep = zeros(B,nchain,rep);

% GRUBIN: chain variance
betamat_chainvar = zeros(p,K0,d-1,nchain,rep);
beta0_chainvar = zeros(p,d-1,nchain,rep);
Bern_K_chainvar = zeros(K0,B,nchain,rep);
tau_chainvar = zeros(B,nchain,rep);

grubin_betamat = zeros(p,K0,d-1,rep);
grubin_beta0 = zeros(p,d-1,rep);
grubin_Bern_K = zeros(K0,B,rep);
grubin_tau = zeros(B,1,rep);


L = (nrun - burn)/thin;

tic;

parfor g = 1:rep
% for g = 1:1
            
    Y_arr = permute(Yn_big_arr(:,:,:,g), [1 3 2]);

    for ch = 1:nchain
    [Q_mat_arr, beta_mat_arr, beta0_arr, Bern_K_arr, tau_arr, K_star_arr, ...
        Q_active_all, beta_mattil_active, Bern_K_permute2, tau_permute2] ...
        = rev_grubin_csp_pos_randini(Y_arr, K, B, nrun, burn, thin, alpha0);
   
    
    K_star_rep(:, g) = [median(K_star_arr(burn+1:thin:end)); iqr(K_star_arr(burn+1:thin:end))];
    
    % z_beta_vec_rep(:, g) = mean(z_beta_vec_arr(:,burn+1:thin:end), 2);
    
    % store other posterior means
    % Q-matrix, permuted
    Q_mat_pomean_rep(:,:,g) = mean(Q_mat_arr(:,:,burn+1:thin:end), 3);
    
%     % A_mat and z posterior mean; A_mat permuted
%     A_mat_pomean_rep(:,:,g) = mean(A_mat_arr(:,:,burn+1:thin:end), 3);
%     z_mat_pomean_rep(:,:,g) = mean(z_tau_mat_arr(:,:,burn+1:thin:end), 3);
    
%     % permute_vec_rep(g,:) = permute_vec;
%     sig2_beta_rep(:,:,g) = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
    
    
    % -- continuous parameters begin here -- %
    % beta_mat, permuted
    beta_mat_pomean_rep(:,:,:,ch,g) = mean(beta_mat_arr(:,:,:,burn+1:thin:end), 4);
    
    % get beta_tilde
    beta_mat_tilde_arr = zeros(p, K, d, nrun);
    for ii=1:nrun
        beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
    end
    
    % beta_mat posterior mean
    beta_mattil_pomean_rep(:,:,:,ch,g) = mean(beta_mat_tilde_arr(:,:,:,burn+1:thin:end), 4);
    
    % beta0 posterior mean
    beta0_pomean_rep(:,:,ch,g) = mean(beta0_arr(:,:,burn+1:thin:end), 3);
    
    % Bern_K posterior mean
    Bern_K_pomean_rep(:,:,ch,g) = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);

    % tau posterior mean
    tau_pomean_rep(:,ch,g) = mean(tau_arr(:,burn+1:thin:end), 2);
    
    
    % GR: within chain variance
    betamat_chainvar(:,:,:,ch,g) = var(beta_mat_tilde_arr(:,1:K0,1:d-1,burn+1:thin:end), [], 4);
    beta0_chainvar(:,:,ch,g) = var(beta0_arr(:,1:d-1,burn+1:thin:end), [], 3);
    Bern_K_chainvar(:,:,ch,g) = var(Bern_K_arr(1:K0,:,burn+1:thin:end), [], 3);
    tau_chainvar(:,ch,g) = var(tau_arr(:,burn+1:thin:end), [], 2);
    
    
    fprintf('Replicate %d Chain %d completed\n', g, ch);
    end
end

for g = 1:rep
    % -- Gelman-Rubin diagnostic -- %
    betamat_B_var = L * squeeze(var(beta_mattil_pomean_rep(:,1:K0,:,:,g), [], 4));
    beta0_B_var = L * squeeze(var(beta0_pomean_rep(:,:,:,g), [], 3));
    BernK_B_var = L * squeeze(var(Bern_K_pomean_rep(1:K0,:,:,g), [], 3));
    tau_B_var = L * squeeze(var(tau_pomean_rep(:,:,g), [], 2));
    
    betamat_W_var = 1/nchain * squeeze(sum(betamat_chainvar(:,1:K0,:,:,g), 4));
    beta0_W_var = 1/nchain * squeeze(sum(beta0_chainvar(:,:,:,g), 3));
    BernK_W_var = 1/nchain * squeeze(sum(Bern_K_chainvar(1:K0,:,:,g), 3));
    tau_W_var = 1/nchain * squeeze(sum(tau_chainvar(:,:,g), 2));
   
    grubin_betamat(:,:,:,g) = (1-1/L) + 1/L * betamat_B_var(:,:,1:d-1) ./ betamat_W_var;
    grubin_beta0(:,:,g) = (1-1/L) + 1/L * beta0_B_var(:,1:d-1) ./ beta0_W_var;
    grubin_Bern_K(:,:,g) = (1-1/L) + 1/L * BernK_B_var ./ BernK_W_var;
    grubin_tau(:,g) = (1-1/L) + 1/L * tau_B_var ./ tau_W_var;
end

% save('mdata_rev_bp_grubin_new.mat')

mytime = toc;



